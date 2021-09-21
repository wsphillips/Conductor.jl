# Classic Hodgkin Huxley neuron with a "current pulse" stimulus
using Conductor, IfElse, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms
import Conductor: Na, K # shorter aliases for Sodium/Potassium

Vₘ = MembranePotential(-65mV) # set default V₀ = -65mV

nav_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

kdr_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         4, name = :n)]

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named Iₑ = IonCurrent(NonIonic)
electrode_pulse = Iₑ ~ IfElse.ifelse(t > 100.0, IfElse.ifelse(t < 200.0, ustrip(Float64, µA, 400pA), 0.0), 0.0)

@named neuron = CompartmentSystem(Vₘ, channels, reversals;
                                  geometry = Sphere(radius = 20µm),
                                  stimuli = [electrode_pulse])

sim = Simulation(neuron, time = 300ms)

using OrdinaryDiffEq

solution = solve(sim, Rosenbrock23())

using Plots
# Plot at 5kHz sampling
plot(solution; plotdensity=Int(1500), size=(1200,800))

