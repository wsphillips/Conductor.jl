
using Conductor, OrdinaryDiffEq, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, mM
import Conductor: Na, K # shorter aliases for Sodium/Potassium

Vₘ = MembranePotential()
Caᵢ = Concentration(Calcium, 0.03mM)
ICa = MembraneCurrent{Calcium}(aggregate = true)

nav_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0 + exp((Vₘ + 25.5)/-5.29)),
         τₘ = 1.32 - 1.26/(1 + exp((Vₘ + 120.0)/-25.0)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0 + exp((Vₘ + 48.9)/5.18)),
         τₕ = (1.34/(1.0 + exp((Vₘ + 62.9)/-10.0)))*(1.5 + 1.0/(1.0 + exp((Vₘ + 34.9)/3.6))))]
         
cas_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((Vₘ+33.0)/-8.1)),
         τₘ = 2.8 + 14.0/(exp((Vₘ+27.0)/10.0) + exp((Vₘ+70.0)/-13.0)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0+exp((Vₘ+60.0)/6.2)),
         τₕ = 120.0 + 300.0/(exp((Vₘ+55.0)/9.0) + exp((Vₘ+65.0)/-16.0)))]

cat_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0 + exp((Vₘ+27.1)/-7.2)),
         τₘ = 43.4 - 42.6/(1.0 + exp((Vₘ+68.1)/-20.5)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0 + exp((Vₘ+32.1)/5.5)),
         τₕ = 210.0 - 179.6/(1.0 + exp((Vₘ+55.0)/-16.9)))]

ka_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((Vₘ+27.2)/-8.7)),
         τₘ = 23.2 - 20.8/(1.0+exp((Vₘ+32.9)/-15.2)),
         p = 3)
    Gate(SteadyStateTau,
         h∞ = 1.0/(1.0+exp((Vₘ+56.9)/4.9)),
         τₕ = 77.2 - 58.4/(1.0+exp((Vₘ+38.9)/-26.5)))]

kca_kinetics = [
    Gate(SteadyStateTau,
         m∞ = (Caᵢ/(Caᵢ+3.0))/(1.0+exp((Vₘ+28.3)/-12.6)),
         τₘ = 180.6 - 150.2/(1.0+exp((Vₘ+46.0)/-22.7)),
         p = 4)]
       
kdr_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((Vₘ+12.3)/-11.8)),
         τₘ = 14.4 - 12.8/(1.0+exp((Vₘ+28.3)/-19.2)),
         p = 4)]

h_kinetics = [
    Gate(SteadyStateTau,
         m∞ = 1.0/(1.0+exp((Vₘ+75.0)/5.5)),
         τₘ = 2/( exp((Vₘ+169.7)/(-11.6)) + exp((Vₘ- 26.7)/(14.3))))]

@named NaV = IonChannel(Sodium, nav_kinetics, 120mS/cm^2) 
@named CaS = IonChannel(Calcium, cas_kinetics, 0mS/cm^2)
@named CaT = IonChannel(Calcium, cat_kinetics, 0mS/cm^2)
@named KA  = IonChannel(Potassium, ka_kinetics, 0mS/cm^2)
@named KCa = IonChannel(Potassium, kca_kinetics, 0mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, 36mS/cm^2)
@named H = IonChannel(Cation, h_kinetics, 10mS/cm^2)
@named leak = PassiveChannel(Leak, 0.3mS/cm^2)


gradients = Equilibria(Pair[Na      =>  50.0mV,
                            K       => -77.0mV,
                            Cation  => (-20mV, :H),
                            Leak    => -54.4mV,
                            Calcium => (500.0)*(8.6174e-5)*(283.15)*(log(max((3000.0/Caᵢ), 0.001)))])

calcium_conversion = AuxConversion(@parameters(τCa = 200.0, Ca∞ = 0.05, f = 0.094),
                                   [D(Caᵢ) ~ (1/τCa)*(-Caᵢ + Ca∞ + f*ICa)])

# problem arises when we have two calcium currents
@named neuron = Soma([NaV, CaS, CaT, KA, KCa, Kdr, H, leak], gradients, aux = [calcium_conversion])

t = 250
sim = Simulation(neuron, time = t*ms)

solution = solve(sim, ImplicitEuler())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(t*5))

