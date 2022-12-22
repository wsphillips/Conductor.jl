
using LinearAlgebra, OrdinaryDiffEq, ModelingToolkit
LinearAlgebra.BLAS.set_num_threads(6)

include(joinpath(@__DIR__, "traub_kinetics.jl"))
include(joinpath(@__DIR__, "pinsky_setup.jl"))

# Base dynamics
@named I_s_holding = Bias(-0.5µA / 0.5)
@named I_d_holding = Bias(0.0µA / (1 - 0.5))

soma_dynamics = HodgkinHuxley([NaV(30mS / cm^2),
                               Kdr(15mS / cm^2),
                               leak(0.1mS / cm^2)],
                               reversals[1:3]);

@named soma = Compartment(Vₘ, soma_dynamics;
                          geometry = Unitless(0.5), # FIXME: should take p param
                          capacitance = capacitance,
                          stimuli = [I_s_holding])

dendrite_dynamics = HodgkinHuxley([KAHP(0.8mS / cm^2),
                                   CaS(10mS / cm^2),
                                   KCa(15mS / cm^2),
                                   leak(0.1mS / cm^2)],
                                  reversals[2:4]);

@named dendrite = Compartment(Vₘ, dendrite_dynamics;
                              extensions = [calcium_conversion],
                              geometry = Unitless(0.5), # FIXME: should take p param
                              capacitance = capacitance,
                              stimuli = [I_d_holding])

@named gc_soma = AxialConductance([Gate(inv(p), name = :ps)], max_g = gc_val)
@named gc_dendrite = AxialConductance([Gate(inv(1 - p), name = :pd)], max_g = gc_val)
topology = MultiCompartmentTopology([soma, dendrite]);
add_junction!(topology, soma, dendrite, (gc_soma, gc_dendrite))
@named mcneuron = MultiCompartment(topology)

###########################################################################################
# Synapse models
###########################################################################################
import Conductor: NMDA, AMPA

@parameters τNMDA = 150.0 τAMPA = 2.0
@variables S(t) = 0.0
@named NMDAChan = SynapticChannel(ConstantValueEvent(S; saturation = 125.0), NMDA,
                                  [Gate(inv(1 + 0.28 * exp(-0.062(Vₘ - 60.0))); name = :e),
                                   Gate(S, [D(S) ~ -S/τNMDA]),
                                   Gate(inv(1 - p), name = :pnmda)],
                                  max_s = 0.027mS)

@named AMPAChan = SynapticChannel(ConstantValueEvent(S; threshold = 20mV), AMPA,
                                  [Gate(S, [D(S) ~ -S/τAMPA]),
                                   Gate(inv(1 - p), name = :pampa)],
                                  max_s = 0.008mS)

ENMDA = EquilibriumPotential(NMDA, 60mV, name = :NMDA)
EAMPA = EquilibriumPotential(AMPA, 60mV, name = :AMPA)
revmap = Dict([NMDAChan => ENMDA, AMPAChan => EAMPA])

# A single neuron was "briefly" stimulated to trigger the network
@named I_pulse = PulseTrain(amplitude = 5.0µA / 0.5,
                            duration = 50ms,
                            delay = 150.0ms,
                            offset = -0.5µA / 0.5)

soma_stimulated = Compartment(Vₘ, deepcopy(soma_dynamics);
                              geometry = Unitless(0.5), # FIXME: should take p param
                              capacitance = capacitance,
                              stimuli = [I_pulse],
                              name = :soma)
mcstim_topology = MultiCompartmentTopology([soma_stimulated, dendrite]);
add_junction!(mcstim_topology, soma_stimulated, dendrite, (gc_soma, gc_dendrite))
@named mcneuron_stim = MultiCompartment(mcstim_topology)

# Need to introduce 10% gca variance as per Pinsky/Rinzel
n_neurons = 100
neuronpopulation = [Conductor.replicate(mcneuron) for _ in 1:n_neurons];
neuronpopulation[4] = mcneuron_stim
topology = NetworkTopology(neuronpopulation, [NMDAChan, AMPAChan]);

using Graphs
nmda_g = random_regular_digraph(n_neurons, fld(n_neurons, 5), dir = :in)
ampa_g = random_regular_digraph(n_neurons, fld(n_neurons, 5), dir = :in)

# We could allow users to supply a lambda/function to map in order to get this behavior
for (i, e) in enumerate(edges(nmda_g))
    add_synapse!(topology, neuronpopulation[src(e)].soma, neuronpopulation[dst(e)].dendrite,
                 NMDAChan, 1.0)
end

for (i, e) in enumerate(edges(ampa_g))
    add_synapse!(topology, neuronpopulation[src(e)].soma, neuronpopulation[dst(e)].dendrite,
                 AMPAChan, 1.0)
end

@named net = NeuronalNetworkSystem(topology, revmap);
t_total = 2000.0
simp = Simulation(net, t_total*ms, return_system = true)
prob = Simulation(net, t_total*ms)
sol = solve(prob, RK4(); adaptive = false, dt = 0.05);

# Pinsky and Rinzel displayed their results as a plot of N neurons over 20mV
indexof(sym, syms) = findfirst(isequal(sym), syms)
dvs = states(simp)
interpolated = sol(1:0.2:t_total,
                   idxs = [indexof(x.soma.Vₘ, dvs) for x in neuronpopulation[5:end]])
abovethold = reduce(hcat, interpolated.u) .> 10.0
using Statistics
final = mean(abovethold, dims = 1)'
using Plots
plot(final)

