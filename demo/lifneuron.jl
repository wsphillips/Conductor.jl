
using Conductor, ModelingToolkit, OrdinaryDiffEq,Graphs, SparseArrays
import Unitful.ms

stim_dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0,
                              threshold = -55.0, resistance = 0.1, stimulus = 500)
@named stim_neuron = CompartmentSystem(stim_dynamics)
sim_neuron = Simulation(stim_neuron, time = 300ms)
# sol_neuron = solve(sim_neuron, Euler(); dt = 0.25);
# lines(sol_neuron[stim_neuron.V] + sol_neuron[stim_neuron.S]*(70))

dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0,
                         resistance = 0.1, stimulus = 0)

neuron_pop = [CompartmentSystem(dynamics) for _ in 1:100]
neuron_pop[5] = stim_neuron
topology = NetworkTopology(neuron_pop, random_regular_digraph(100, 20, dir=:in);
                           default_weight=90.0);

@named network = NeuronalNetworkSystem(topology)
sim_network = Simulation(network, time = 300ms);
@timed sol_network = solve(sim_network, Euler(), adaptive = false, dt = 1);

# kills the compiler for large neuron populations, so create it manually for now
deduplicated_grid = sol_network(0.0:1.0:300.0)
# with a passthrough state, no problem
series = [deduplicated_grid[network[x].V]' for x in 1:50]
# but with an observable (RGF eval for each subsystem) can cause lockup or total crash
deduplicated_grid[network[1].S]
series = [deduplicated_grid[network[x].S] for x in 1:100]

spike_map = sparse(hcat((sol_network[network[x].S] for x in 1:50)...))

spike_map = sparse(map(sol_network[1:2:200,:]) do x
             x >= -55.0
         end)

using GLMakie, Makie
spy(spike_map, markersize = 4, marker = :circle)
# plot(sol_network, vars = [network[x].S for x in 1:100])
# plot(sol_network[network[5].V])
# 

