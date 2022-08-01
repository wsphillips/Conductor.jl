
using Conductor, ModelingToolkit, OrdinaryDiffEq, Plots
import Unitful.ms

dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0, resistance = 0.1, stimulus = 500)
    @named neuron = CompartmentSystem(dynamics)
    neuron_pop = [CompartmentSystem(dynamics) for _ in 1:10]
    topology = NetworkTopology(neuron_pop);

    topology[3,4] = 1.0
    topology[4,6] = 0.5

@named network = NeuronalNetworkSystem(topology)

sim = Simulation(neuron, time = 300ms, return_system = false)
@time sol = solve(sim, Rosenbrock23(), dt = 0.5);
plot(sol[neuron.V] + sol[neuron.S]*(70)) # Superimpose spikes on top of voltage
plot(sol)

sim2 = Simulation(network, time = 300ms, return_system = false)
@time sol = solve(sim2, Rosenbrock23(), dt = 0.5);
plot(sol)

