# Classic Hodgkin Huxley neuron with a "current pulse" stimulus
using Conductor, Unitful, ModelingToolkit, OrdinaryDiffEq, Plots
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms
import Conductor: Na, K # shorter aliases for Sodium/Potassium

Vₘ = MembranePotential(-65mV)

nav_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)];

kdr_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4, name = :n)];

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named Iₑ = PulseTrain(amplitude = 400.0pA, duration = 100ms, delay = 100ms)

dynamics = HodgkinHuxley(Vₘ, channels, reversals;
                         geometry = Sphere(radius = 20µm),
                         stimuli = [Iₑ]);

@named neuron = Compartment(dynamics)
sim = Simulation(neuron, time = 300ms)
solution = solve(sim, Rosenbrock23(), abstol = 0.01, reltol = 0.01, dtmax = 100.0);
plot(solution; size=(1200,800))

