
using Conductor, ModelingToolkit, OrdinaryDiffEq,Graphs, SparseArrays, Plots
import Unitful.ms

stim_dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0,
                              threshold = -55.0, resistance = 0.1, stimulus = 200.01)
@named stim_neuron = CompartmentSystem(stim_dynamics)
sim_neuron = Simulation(stim_neuron, time = 1000ms)
sol_neuron = solve(sim_neuron, Euler(); dt = 0.25);
plot(sol_neuron[stim_neuron.V] + sol_neuron[stim_neuron.S]*(70)) # single neuron ok


dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0,
                         resistance = 0.1, stimulus = 0)
neuron_pop = [CompartmentSystem(dynamics) for _ in 1:10]
neuron_pop[5] = stim_neuron
topology = NetworkTopology(neuron_pop, random_regular_digraph(10, 4, dir=:in);
                           default_weight=100.0);

# NOTE: This is configured to causes network-wide spiking but it will crash even with
# completely silent network (set `default_weight = 80.0` for subthreshold connectivity)
@named network = NeuronalNetworkSystem(topology)
sim_network = Simulation(network, time = 300ms, return_system=true);

# I have tried several solvers and callback triggering mechanisms. Euler is currently least
# likely to crash out.
@timed sol_network = solve(sim_network, Euler(), dt=1); # careful--this might end your day

deduplicated_grid = sol_network(0.0:1.0:300.0, idxs=[network[x].S for x in 1:100]

# Symbolic solution indexing with observables will also cause issues.
# It hangs or crashes for large neuron populations, so handle data manually for now
# Example, with a passthrough state, no problem
#series = [deduplicated_grid[network[x].V]' for x in 1:50]
# ...but with an observable (RGF eval for each subsystem) can cause lockup or total crash
#deduplicated_grid[network[1].S] # just one usually works, but be careful
#series = [deduplicated_grid[network[x].S] for x in 1:100] # will 100% hang or terminate process
#spike_map = sparse(hcat((sol_network[network[x].S] for x in 1:50)...))

using GLMakie, Makie
spy(spike_map, markersize = 4, marker = :circle)
# plot(sol_network, vars = [network[x].S for x in 1:100])
# plot(sol_network[network[5].V])
# 

