
using Conductor, ModelingToolkit, OrdinaryDiffEq, Plots, Graphs
import Unitful.ms

stim_dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0, resistance = 0.1, stimulus = 500)
@named stim_neuron = CompartmentSystem(stim_dynamics)

sim = Simulation(stim_neuron, time = 300ms, return_system = false)
@time sol = solve(sim, Euler(), dt = 0.5, adaptive = false);
plot(sol[stim_neuron.V] + sol[stim_neuron.S]*(70)) # Superimpose spikes on top of voltage

dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0, resistance = 0.1, stimulus = 0)
neuron_pop = [CompartmentSystem(dynamics) for _ in 1:10]
neuron_pop[5] = stim_neuron
topology = NetworkTopology(neuron_pop, random_regular_digraph(10, 4, dir=:in); default_weight=5000.0);

@named network = NeuronalNetworkSystem(topology)

sim2 = Simulation(network, time = 300ms, return_system = false);
@time sol2 = solve(sim2, Euler(), dt = 0.5, adaptive=false);
plot(sol2, vars = [network[x].S for x in 1:10])

