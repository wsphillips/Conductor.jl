
using Conductor, ModelingToolkit
import Unitful.ms
dynamics = Conductor.LIF(-75.0, tau_membrane = 10.0, tau_synaptic = 10.0, threshold = -55.0, resistance = 0.1, stimulus = 500)

@named neuron = CompartmentSystem(dynamics)

sim = Simulation(neuron, time = 300ms, return_system = false)

using OrdinaryDiffEq, Plots

sol = solve(sim, Rosenbrock23())
plot(sol)
