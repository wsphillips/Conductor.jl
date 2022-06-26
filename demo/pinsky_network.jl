
cd("Conductor.jl/demo")

include(joinpath(@__DIR__, "traub_kinetics.jl"))

import Unitful: µF, pA, µA, nA, µS

@parameters ϕ = 0.13 β = 0.075 p = 0.5

@named calcium_conversion = ODESystem([D(Caᵢ) ~ -ϕ*ICa - β*Caᵢ]);

reversals = Equilibria(Pair[Sodium    =>  120.0mV,
                            Potassium =>  -15.0mV,
                            Leak      =>    0.0mV,
                            Calcium   =>  140.0mV]);

capacitance = 3.0µF/cm^2
gc_val = 2.1mS/cm^2

# Pinsky modifies NaV to have instantaneous activation, so we can ignore tau
pinsky_nav_kinetics = [convert(Gate{SteadyState}, nav_kinetics[1]), nav_kinetics[2]]
@named NaV = IonChannel(Sodium, pinsky_nav_kinetics) 

# No inactivation term for calcium current in Pinsky model
pinsky_ca_kinetics = [ca_kinetics[1]]
@named CaS = IonChannel(Calcium, pinsky_ca_kinetics)

is_val = ustrip(Float64, µA, 0.75µA)/p
@named Iₛ = IonCurrent(NonIonic, is_val, dynamic = false)
soma_holding = Iₛ ~ is_val

id_val = ustrip(Float64, µA, -0.5µA)/(1-p)
@named I_d = IonCurrent(NonIonic, id_val, dynamic = false)
dendrite_holding = I_d ~ id_val

@named soma = Compartment(Vₘ,
                         [NaV(30mS/cm^2),
                          Kdr(15mS/cm^2),
                          leak(0.1mS/cm^2)],
                          reversals[1:3],
                          geometry = Unitless(0.5), # FIXME: should take p param
                          capacitance = capacitance,
                          stimuli = [soma_holding])

@named dendrite = Compartment(Vₘ,
                             [KAHP(0.8mS/cm^2),
                              CaS(10mS/cm^2),
                              KCa(15mS/cm^2),
                              leak(0.1mS/cm^2)],
                             reversals[2:4],
                             geometry = Unitless(0.5), # FIXME: should take p param
                             capacitance = capacitance,
                             extensions = [calcium_conversion])

@named gc_soma = AxialConductance([Gate(SteadyState, inv(p), name = :ps)], max_g = gc_val)
@named gc_dendrite = AxialConductance([Gate(SteadyState, inv(1-p), name = :pd)], max_g = gc_val)
soma2dendrite = Junction(soma => dendrite, gc_soma, symmetric = false);
dendrite2soma = Junction(dendrite => soma, gc_dendrite, symmetric = false);
@named mcneuron = MultiCompartment([soma2dendrite, dendrite2soma])

###########################################################################################
# Synapse models
###########################################################################################
import Conductor: NMDA, AMPA, HeavisideSum

@named NMDAChan = SynapticChannel(NMDA,
                [Gate(SteadyState, inv(1 + 0.28*exp(-0.062(Vₘ - 60.))); name = :e),
                 Gate(HeavisideSum, threshold = 10mV, saturation = 150; name = :S),
                 Gate(SteadyState, inv(1-p), name = :pnmda)],
                 max_s = 0.014mS, aggregate = true)

@named AMPAChan = SynapticChannel(AMPA,
                                [Gate(HeavisideSum, threshold = 20mV, saturation = 2;
                                      name = :u),
                                 Gate(SteadyState, inv(1-p), name = :pampa)],
                                 max_s = 0.0045mS, aggregate = true)

ENMDA = EquilibriumPotential(NMDA, 60mV, name = :NMDA)
EAMPA = EquilibriumPotential(AMPA, 60mV, name = :AMPA)

@time neuronpopulation = [Conductor.replicate(mcneuron) for n in 1:100];

allsynapses = Set{Synapse}()

# This will go much faster using something like SimpleDiGraph--currently way too much
# allocation
@time begin
for neuron in neuronpopulation

    ampa_synapses = Set{Synapse}()
    nmda_synapses = Set{Synapse}()
    while length(nmda_synapses) < 20  
        syn1 = Synapse(rand(neuronpopulation).soma => neuron.dendrite, NMDAChan, ENMDA)
        push!(nmda_synapses, syn1)
    end 

    while length(ampa_synapses) < 20
        syn2 = Synapse(rand(neuronpopulation).soma => neuron.dendrite, AMPAChan, EAMPA)
        push!(ampa_synapses, syn2)
    end
    union!(allsynapses, ampa_synapses, nmda_synapses);
end
end
allsynapses = collect(allsynapses);

# this is only fast because the cost is amortized
@time net = NeuronalNetworkSystem(allsynapses);

# the steps include a slow conversion: convert(ODESystem, net)
# structural_simplify is only a minor cost (~12 seconds for 100 neurons with 4000 synapses)
prob = Simulation(net, time = 5000ms)
using OrdinaryDiffEq
# on the first run, over 40% of the time is single thread compilation time. Without
# compilation, the solution is about 200 seconds for 5 seconds of tspan
sol = solve(prob, RadauIIA5())
using Plots
plot(sol)

