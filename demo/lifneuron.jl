
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

neuron_pop = [CompartmentSystem(dynamics) for _ in 1:50]
neuron_pop[5] = stim_neuron
topology = NetworkTopology(neuron_pop, random_regular_digraph(50, 20, dir=:in);
                           default_weight=110.0);

@named network = NeuronalNetworkSystem(topology)
sim_network = Simulation(network, time = 300ms);
@timed sol_network = solve(sim_network, Euler(), adaptive = false, dt = 1);

# kills the compiler for large neuron populations, so create it manually for now
# spike_map = sparse(vcat((sol_network[network[x].S]' for x in 1:50)...))
 spike_map = sparse(map(sol_network[1:2:100,:]) do x
             x >= -55.0
         end)


using GLMakie, Makie
spy(spike_map, markersize = 4, marker = :circle)
# plot(sol_network, vars = [network[x].S for x in 1:100])
# plot(sol_network[network[5].V])
# 

