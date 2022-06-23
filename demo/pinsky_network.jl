
cd("/home/wsphil/git/Conductor.jl/demo")
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

neuronpopulation = [Conductor.replicate(mcneuron) for n in 1:10];

allsynapses = Synapse[]

for neuron in neurons
    for _ in 1:2
        syn1 = Synapse(rand(neuronpopulation).soma => neuron.dendrite, NMDAChan, ENMDA)
        syn2 = Synapse(rand(neuronpopulation).soma => neuron.dendrite, AMPAChan, EAMPA)
        push!(allsynapses, syn1)
        push!(allsynapses, syn2)
    end
end

net = NeuronalNetworkSystem(allsynapses)

