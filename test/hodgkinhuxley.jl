
using Conductor, IfElse, OrdinaryDiffEq, Plots, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms
import Conductor: Na, K # shorter aliases for Sodium/Potassium

Vₘ = MembranePotential()

nav_kinetics = [
    Gate(AlphaBetaRates,
         αₘ = IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         βₘ = 4.0*exp(-(Vₘ + 65.0)/18.0),
         p = 3)
    Gate(AlphaBetaRates,
         αₕ = 0.07*exp(-(Vₘ+65.0)/20.0),
         βₕ = 1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)))]

kdr_kinetics = [
    Gate(AlphaBetaRates,
         αₙ = IfElse.ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         βₙ = 0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4)]

@named NaV = IonChannel(Sodium, nav_kinetics, 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, 36mS/cm^2)
@named leak = PassiveChannel(Leak, 0.3mS/cm^2)

# Equilibrium potentials are a implicit description of a ion concentration gradient
gradients = Equilibria([Na   =>  50.0mV,
                        K    => -77.0mV,
                        Leak => -54.4mV])

area = 4*pi*(20µm)^2
pulse(t, current) = 100. < t < 200. ? ustrip(Float64, µA, 400pA) : 0.0
@register pulse(a,b)

@named neuron = Soma([NaV,Kdr,leak], gradients, stimulus = pulse, area = ustrip(Float64, cm^2, area));

t = 300 
sim = Simulation(neuron, time = t*ms)

solution = solve(sim, Rosenbrock23())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(t*5))

