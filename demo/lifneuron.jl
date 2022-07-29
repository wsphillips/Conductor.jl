
using Conductor, ModelingToolkit
import Unitful.ms
dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0, resistance = 0.1, stimulus = 500)

@named neuron = CompartmentSystem(dynamics)

neuron_pop = [CompartmentSystem(dynamics) for _ in 1:10]
net = NetworkTopology(neuron_pop);

sim = Simulation(neuron, time = 300ms, return_system = false)

using OrdinaryDiffEq, Plots

# LIF neuron -- (fixed time step only for now)
sol = solve(sim, Rosenbrock23(), dt = 0.5);

# Superimpose spikes on top of voltage
plot(sol[neuron.V] + sol[neuron.S]*(70))


