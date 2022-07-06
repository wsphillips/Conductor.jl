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
#pinsky_nav_kinetics = [convert(Gate{SteadyState}, nav_kinetics[1]), nav_kinetics[2]]
#@named NaV = IonChannel(Sodium, pinsky_nav_kinetics) 

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
                             extensions = [calcium_conversion],
                             stimuli = [dendrite_holding])


@named gc_soma = AxialConductance([Gate(SimpleGate, inv(p), name = :ps)],
                                  max_g = gc_val)
@named gc_dendrite = AxialConductance([Gate(SimpleGate, inv(1-p), name = :pd)],
                                      max_g = gc_val)

topology = Conductor.MultiCompartmentTopology([soma, dendrite]);

Conductor.add_junction!(topology, soma,  dendrite, gc_soma, symmetric = false)
Conductor.add_junction!(topology, dendrite,  soma, gc_dendrite, symmetric = false)

@named mcneuron = MultiCompartment(topology)

###########################################################################################
# Synapse models
###########################################################################################
import Conductor: NMDA, AMPA, HeavisideSum

@named NMDAChan = SynapticChannel(NMDA,
                [Gate(SimpleGate, inv(1 + 0.28*exp(-0.062(Vₘ - 60.))); name = :e),
                 Gate(HeavisideSum, threshold = 10mV, saturation = 150; name = :S),
                 Gate(SimpleGate, inv(1-p), name = :pnmda)],
                 max_s = 0.014mS, aggregate = true)

@named AMPAChan = SynapticChannel(AMPA,
                                [Gate(HeavisideSum, threshold = 20mV, saturation = 2;
                                      name = :u),
                                 Gate(SimpleGate, inv(1-p), name = :pampa)],
                                 max_s = 0.0045mS, aggregate = true)

ENMDA = EquilibriumPotential(NMDA, 60mV, name = :NMDA)
EAMPA = EquilibriumPotential(AMPA, 60mV, name = :AMPA)

# Need to introduce 10% gca variance as per Pinsky/Rinzel
neuronpopulation = Vector{MultiCompartmentSystem}(undef, 100)

# introduce calcium channel conductance variation
for i in eachindex(neuronpopulation)
    @named dendrite = Compartment(Vₘ,
                                 [KAHP(0.8mS/cm^2),
                                  CaS(rand(9.0:0.1:10.0)mS/cm^2),
                                  KCa(15mS/cm^2),
                                  leak(0.1mS/cm^2)],
                                 reversals[2:4],
                                 geometry = Unitless(0.5),
                                 capacitance = capacitance,
                                 extensions = [calcium_conversion])
    
    soma2dendrite = deepcopy(Junction(soma => dendrite, gc_soma, symmetric = false));
    dendrite2soma = deepcopy(Junction(dendrite => soma, gc_dendrite, symmetric = false));
    neuronpopulation[i] = MultiCompartment([soma2dendrite, dendrite2soma]; name = Conductor.namegen(:mcneuron))
   
end

using Graphs
allsynapses = Vector{Synapse}(undef, 4000)

# Generating topology with Graphs.jl is trivial.
# The getproperty calls for soma/dendrite are far too expensive here.
@time begin
nmda_g = random_regular_digraph(100, 20, dir=:in)
ampa_g = random_regular_digraph(100, 20, dir=:in)
for (i, e) in enumerate(edges(nmda_g))
    allsynapses[i] = Synapse(neuronpopulation[src(e)].soma => neuronpopulation[dst(e)].dendrite, NMDAChan, ENMDA)
end
for (i, e) in enumerate(edges(ampa_g))
    allsynapses[2000 + i] = Synapse(neuronpopulation[src(e)].soma => neuronpopulation[dst(e)].dendrite, AMPAChan, EAMPA)
end
end

# this is only fast because the cost is amortized
@time net = NeuronalNetworkSystem(allsynapses);

# the steps include a way too slow conversion: convert(ODESystem, net)
# relatively, call to structural_simplify is cheap

simp = Simulation(net, time = 5000ms, return_system=true)
prob = ODAEProblem(simp, [], (0,5000.));
# on the first run, over 40% of the time is single threaded compilation time. Without
# compilation, solving takes roughly 200 seconds (6-8 cores) for 5 seconds of tspan.

@time sol = solve(prob, RadauIIA5());
plot(sol)

