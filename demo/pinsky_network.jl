

using LinearAlgebra, OrdinaryDiffEq
LinearAlgebra.BLAS.set_num_threads(6)

include(joinpath(@__DIR__, "traub_kinetics.jl"))
include(joinpath(@__DIR__, "pinsky_setup.jl"))

# Base dynamics
@named I_s_holding = IonCurrent(NonIonic, -0.5µA, dynamic = false)
@named Iₛ = IonCurrent(NonIonic, -0.5µA, dynamic = false)
soma_holding = Iₛ ~ I_s_holding/p

@named I_d_holding = IonCurrent(NonIonic, 0.0µA, dynamic = false)
@named I_d = IonCurrent(NonIonic, 0.0µA, dynamic = false)
dendrite_holding = I_d ~ I_d_holding/(1-p)

soma_dynamics = HodgkinHuxley(Vₘ,
                         [NaV(30mS/cm^2),
                          Kdr(15mS/cm^2),
                          leak(0.1mS/cm^2)],
                          reversals[1:3];
                          geometry = Unitless(0.5), # FIXME: should take p param
                          capacitance = capacitance,
                          stimuli = [soma_holding])

@named soma = Compartment(soma_dynamics)

dendrite_dynamics = HodgkinHuxley(Vₘ,
                             [KAHP(0.8mS/cm^2),
                              CaS(10mS/cm^2),
                              KCa(15mS/cm^2),
                              leak(0.1mS/cm^2)],
                             reversals[2:4],
                             geometry = Unitless(0.5), # FIXME: should take p param
                             capacitance = capacitance,
                             stimuli = [dendrite_holding])

@named dendrite = Compartment(dendrite_dynamics,
                              extensions = [calcium_conversion])

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
                 Gate(HeavisideSum; threshold = 10mV, decay = 150, saturation = 125, name = :S),
                 Gate(SimpleGate, inv(1-p), name = :pnmda)],
                 max_s = 0.028mS, aggregate = true)

@named AMPAChan = SynapticChannel(AMPA,
                                [Gate(HeavisideSum, threshold = 20mV, decay = 2;
                                      name = :u),
                                 Gate(SimpleGate, inv(1-p), name = :pampa)],
                                 max_s = 0.018mS, aggregate = true)

ENMDA = EquilibriumPotential(NMDA, 60mV, name = :NMDA)
EAMPA = EquilibriumPotential(AMPA, 60mV, name = :AMPA)
revmap = Dict([NMDAChan => ENMDA, AMPAChan => EAMPA])

# A single neuron was "briefly" stimulated to trigger the network
@named I_pulse = IonCurrent(NonIonic, 5.0µA; dynamic = false)
@named I_holding = IonCurrent(NonIonic, -0.5µA; dynamic = false)
@parameters t_on = 150.0 [unit=ms] t_off = 175.0 [unit=ms]
@named Istim = IonCurrent(NonIonic, -0.5µA)
soma_stim = Istim ~ ifelse((t > t_on) & (t < t_off), I_pulse/p, I_holding/p)

soma_stim_dynamics = HodgkinHuxley(Vₘ,
                                   [NaV(30mS/cm^2),
                                    Kdr(15mS/cm^2),
                                    leak(0.1mS/cm^2)],
                                    reversals[1:3];
                                    geometry = Unitless(0.5), # FIXME: should take p param
                                    capacitance = capacitance,
                                    stimuli = [soma_stim])

soma_stimulated = Compartment(soma_stim_dynamics; name=:soma)
mcstim_topology = Conductor.MultiCompartmentTopology([soma_stimulated, dendrite]);
Conductor.add_junction!(mcstim_topology, soma_stimulated,  dendrite, gc_soma, symmetric = false)
Conductor.add_junction!(mcstim_topology, dendrite,  soma_stimulated, gc_dendrite, symmetric = false)
@named mcneuron_stim = MultiCompartment(mcstim_topology)

# Need to introduce 10% gca variance as per Pinsky/Rinzel
neuronpopulation = [Conductor.replicate(mcneuron) for _ in 1:100];
neuronpopulation[4] = mcneuron_stim
topology = NetworkTopology(neuronpopulation, [NMDAChan, AMPAChan]);

using Graphs
nmda_g = random_regular_digraph(100, 20, dir=:in)
ampa_g = random_regular_digraph(100, 20, dir=:in)

# We could allow users to supply a lambda/function to map in order to get this behavior
for (i, e) in enumerate(edges(nmda_g))
    add_synapse!(topology, neuronpopulation[src(e)].soma, neuronpopulation[dst(e)].dendrite, NMDAChan)
end

for (i, e) in enumerate(edges(ampa_g))
    add_synapse!(topology, neuronpopulation[src(e)].soma, neuronpopulation[dst(e)].dendrite, AMPAChan)
end

@named net = NeuronalNetworkSystem(topology, revmap);
simp = Simulation(net, time = 2000.0ms, return_system = true)
prob = Simulation(net, time = 2000.0ms)
@time sol = solve(prob, RK4());

# Pinsky and Rinzel displayed their results as a plot of N neurons over 20mV
indexof(sym,syms) = findfirst(isequal(sym),syms)
dvs = states(simp)
interpolated = sol(1:0.2:2000.0, idxs=[indexof(x.soma.Vₘ, dvs) for x in neuronpopulation[5:end]])
abovethold = reduce(hcat, interpolated.u) .> 10.0
using Statistics
final = mean(abovethold, dims=1)'
using Plots
plot(final) # looks correct but only for less than 1 second of simulation time.
