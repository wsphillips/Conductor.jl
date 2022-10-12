
using Conductor, ModelingToolkit, OrdinaryDiffEq, Graphs, SparseArrays, Plots
import Unitful.ms

stim_dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0,
                              threshold = -55.0, resistance = 0.1, stimulus = 220.01)
@named stim_neuron = CompartmentSystem(stim_dynamics)

sim_neuron = Simulation(stim_neuron, time = 1000ms)
sol_neuron = solve(sim_neuron, Euler(); dt = 0.25, saveat=0.25,
                   save_idxs=1);
data = Array(sol_neuron)[1:4000]

function loss(p)
    _prob = remake(sim_neuron, p=p);
    _sol = solve(_prob, Euler(), dt = 0.25, saveat=0.25, save_idxs=1, sensealg = ForwardDiffSensitivity());

    

    l = sum(abs2, Array(_sol)[1:4000] .- data);
    return l
end

using SciMLSensitivity, Optimization, OptimizationOptimJL, ForwardDiff
using OptimizationMultistartOptimization
p0 = [0.0,0.0,-80.0,0.0,0.0,-50.0]

#ReverseDiff.gradient(loss, p0)
#ForwardDiff.gradient(loss, p0)
#prob = remake(sim_neuron, p=p0)
#S = observed(stim_neuron)[1].lhs # hack until bug-fix in MTK
#Plots.plot(sol_neuron[stim_neuron.V] + sol_neuron[S]*(70))

adtype = Optimization.AutoForwardDiff()
optf = Optimization.OptimizationFunction((x,p)->loss(x), adtype)
prob = Optimization.OptimizationProblem(optf, p0, lb = [0.0,0.0,-100.0, 0.0, 0.0, -100.0],
                                                  ub = [0.01,1000.0,100.0,1000.0,1000.0,100.0])

callback = function (p,l,pred) #callback function to observe training
    display(l)
    #plot(solve(remake(sim_neuron, p=p), Euler(), dt=0.25)) |> display
    return false # Tell it to not halt the optimization. If return true, then optimization stops
end

sol = solve(prob, MultistartOptimization.TikTak(1000), LBFGS(), maxiters=50000)
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

