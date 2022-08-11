
using Conductor, ModelingToolkit, OrdinaryDiffEq,Graphs, SparseArrays, Plots
import Unitful.ms

stim_dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0,
                              threshold = -55.0, resistance = 0.1, stimulus = 220.01)
@named stim_neuron = CompartmentSystem(stim_dynamics)

sim_neuron = Simulation(stim_neuron, time = 1000ms)
sol_neuron = solve(sim_neuron, Euler(); dt = 0.25);

S = observed(stim_neuron)[1].lhs # hack until bug-fix in MTK
Plots.plot(sol_neuron[stim_neuron.V] + sol_neuron[S]*(70))

# Network simulation, constant spiking
n = 100 # number of neurons
w = 400.0 # synaptic weight (homogenous)
t_total = 500.0 # total time of simulation
dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0,
                         resistance = 0.1, stimulus = 0);
neuron_pop = [CompartmentSystem(dynamics) for _ in 1:n]; # should use Population instead
neuron_pop[5] = stim_neuron;
topology = NetworkTopology(neuron_pop, random_regular_digraph(n, 20, dir=:in);
                           default_weight=w);
@named network = NeuronalNetworkSystem(topology);
sim_network = Simulation(network, time = t_total*ms);
@time sol_network = solve(sim_network, Rosenbrock23(), tstops=0.0:1.0:t_total);
dedup_grid = sol_network(0.0:1.0:300.0) # deduplicate timestep saving from callbacks
# NOTE: Symbolic solution indexing with can cause issues, so handle data manually for now
spike_map = sparse(map(x -> x >= -55.0, dedup_grid[1:2:2*n,:]))

# To see a raster plot:
using GLMakie, Makie
Makie.spy(spike_map, markersize = 4, marker = :circle)

