
using Conductor, IfElse, OrdinaryDiffEq, Plots, Unitful
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

@named neuron = Soma([NaV,Kdr,leak], gradients, applied = 225nA, radius = 20µm)

t = 250
sim = Simulation(neuron, time = t*ms)

solution = solve(sim, Rodas5())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(t*5))

