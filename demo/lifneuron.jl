
using Conductor, ModelingToolkit, OrdinaryDiffEq, Makie, GLMakie, Graphs, SparseArrays
import Unitful.ms

stim_dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0,
                              threshold = -55.0, resistance = 0.1, stimulus = 500)
@named stim_neuron = CompartmentSystem(stim_dynamics)
sim_neuron = Simulation(stim_neuron, time = 300ms)
@time sol_neuron = solve(sim_neuron, Euler(), dt = 0.25);
plot(sol_neuron[stim_neuron.V] + sol_neuron[stim_neuron.S]*(70))

dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0,
                         resistance = 0.1, stimulus = 0)

neuron_pop = [CompartmentSystem(dynamics) for _ in 1:50]
neuron_pop[5] = stim_neuron
topology = NetworkTopology(neuron_pop, random_regular_digraph(50, 10, dir=:in);
                           default_weight=100.0);

@named network = NeuronalNetworkSystem(topology)
sim_network = Simulation(network, time = 300ms);
@time sol_network = solve(sim_network, ImplicitEuler(), adaptive = false, dt = 0.5, saveat=0.5);

even_grid = sol_network(0.0:0.5:300.0)

spike_map = sparse(vcat( (even_grid[network[x].S]' for x in 1:50)...))
spy(spike_map, markersize = 4, marker = :circle)
plot(sol_network, vars = [network[x].S for x in 1:100])
plot(sol_network[network[5].V])

#simp_network = Simulation(network, time = 300ms, return_system = true);
#pre, sol_states = ModelingToolkit.get_substitutions_and_solved_states(simp_network)
#sym_cb = ModelingToolkit.get_discrete_events(simp_network)
#big_condition = ModelingToolkit.condition(sym_cb[1])
#needed_obseqs = filter(x -> ∉(x.lhs, keys(sol_states.rewrites)), observed(simp_network))
#obs_lhss = [x.lhs for x in needed_obseqs]
#obs_rhss = [x.rhs for x in needed_obseqs]
#cb_subs = Dict(zip(obs_lhss, obs_rhss))
#revised_condition = substitute(big_condition, cb_subs)
#build_function([big_condition], states(simp_network), Conductor.t,
#               parameters(network); expression = Val{true}, states = sol_states,
#               postprocess_fbody = pre)
#compiled_cb = sim2.kwargs[:callback].discrete_callbacks[1]

