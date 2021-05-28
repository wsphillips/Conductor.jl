
using Conductor, OrdinaryDiffEq, Plots, Unitful
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms
import Conductor: Na, K # shorter aliases for Sodium/Potassium

Vₘ = MembranePotential()

nav_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0 + exp((V + 25.5)/-5.29)),
         τₘ = 1.32 - 1.26/(1 + exp((V + 120.0)/-25.0)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0 + exp((V + 48.9)/5.18)),
         τₕ = (1.34/(1.0 + exp((V + 62.9)/-10.0)))*(1.5 + 1.0/(1.0 + exp((V + 34.9)/3.6))))]
         
cas_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((V+33.0)/-8.1)),
         τₘ = 2.8 + 14.0/(exp((V+27.0)/10.0) + exp((V+70.0)/-13.0)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0+exp((V+60.0)/6.2)),
         τₕ = 120.0 + 300.0/(exp((V+55.0)/9.0) + exp((V+65.0)/-16.0)))]

cat_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0 + exp((V+27.1)/-7.2)),
         τₘ = 43.4 - 42.6/(1.0 + exp((V+68.1)/-20.5)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0 + exp((V+32.1)/5.5)),
         τₕ = 210.0 - 179.6/(1.0 + exp((V+55.0)/-16.9)))]

ka_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((V+27.2)/-8.7)),
         τₘ = 23.2 - 20.8/(1.0+exp((V+32.9)/-15.2)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0+exp((V+56.9)/4.9)),
         τₕ = 77.2 - 58.4/(1.0+exp((V+38.9)/-26.5)))]

kca_kinetics = [
    Gate(SteadyStateTau,
         m∞ = (Ca/(Ca+3.0))/(1.0+exp((V+28.3)/-12.6)),
         τₘ = 180.6 - 150.2/(1.0+exp((V+46.0)/-22.7)),
         p = 4)]
       
kdr_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((V+12.3)/-11.8)),
         τₘ = 14.4 - 12.8/(1.0+exp((V+28.3)/-19.2)),
         p = 4)]

h_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((V+75.0)/5.5)),
         τₘ = 2/( exp((V+169.7)/(-11.6)) + exp((V- 26.7)/(14.3))))]

@named NaV = IonChannel(Sodium, nav_kinetics, 120mS/cm^2) 
@named CaS = IonChannel(Calcium, cas_kinetics, 0mS/cm^2)
@named CaT = IonChannel(Calcium, cat_kinetics, 0mS/cm^2)
@named KA  = IonChannel(Potassium, ka_kinetics, 0mS/cm^2)
@named KCa = IonChannel(Potassium, kca_kinetics, 0mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, 0mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, 36mS/cm^2)
@named leak = PassiveChannel(Leak, 0.3mS/cm^2)

# Equilibrium potentials are a implicit description of a ion concentration gradient
gradients = Equilibria([Na   =>  50.0mV,
                        K    => -77.0mV,
                        Leak => -54.4mV])

# ECa = (500.0)*(8.6174e-5)*(283.15)*(log(max((3000.0/Ca), 0.001)))
# f = 0.094; summed_calcium_flux -> sum of all calcium currents
# need params: τCa Ca∞
# D(Ca) ~ (1/τCa)*(-Ca + Ca∞ + f*summed_calcium_flux)
@named neuron = Soma([NaV,Kdr,leak], gradients, applied = 225nA, radius = 20µm)

t = 250
sim = Simulation(neuron, time = t*ms)

solution = solve(sim, Rodas5())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(t*5))

