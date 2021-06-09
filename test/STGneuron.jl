
using Conductor, OrdinaryDiffEq, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, mM, nM, µM
import Conductor: Na, K # shorter aliases for Sodium/Potassium

Vₘ = MembranePotential()
Caᵢ = Concentration(Calcium, 0.05µM)
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

# TODO: Check if calcium concentration should be in mM here...it seems like it.
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

@named NaV = IonChannel(Sodium, nav_kinetics, 100mS/cm^2) 
@named CaT = IonChannel(Calcium, cat_kinetics, 0mS/cm^2)
@named CaS = IonChannel(Calcium, cas_kinetics, 4mS/cm^2)
@named KA  = IonChannel(Potassium, ka_kinetics, 0mS/cm^2)
@named KCa = IonChannel(Potassium, kca_kinetics, 15mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, 50mS/cm^2)
@named H = IonChannel(Cation, h_kinetics, 0.02mS/cm^2)
@named leak = PassiveChannel(Leak, 0.03mS/cm^2)


gradients = Equilibria(Pair[Na      =>  50.0mV,
                            K       => -80.0mV,
                            Cation  => (-20mV, :H),
                            Leak    => -50mV,
                            # x1000 because answer natively in volts (we need mV)
                            # TODO: write a nernst function so we don't have to do this
                            Calcium => 1000*((-8.314*283.15)/(2 * 96485.365))*log(Caᵢ/3000.0)])


# TODO: combine the calcium equilibrium and "conversion" into one clean data structure/MTK system
# internal current units are µA (mS * mV) so we correct by x1000 to match Liu conversion factor
# Perhaps make a "per unit area" version of the conversion factor? It's based on using
# capacitance as an implicit definition of membrane area...when capacitance
# technically could vary independently (although 99% of the time people just assume "1")
# Prinz's units here irk me too. The time constant/f value seem to be shifted by a factor of
# 10, while _not_ scaling internal calcium/steady state calcium. This manipulation
# reproduces the models in her papers but doesn't follow the supposed function described by
# Liu...
# Note also that the 14.96 number comes from dividing 0.94/capacitance * 10 (presumably because
# τCa is 200ms, which is 10x the 20ms used in Liu's paper...
calcium_conversion = AuxConversion(@parameters(τCa = 200.0, Ca∞ = 0.05, f = 14.96),
                                   [D(Caᵢ) ~ (1/τCa)*(-f*(ICa*1000) - Caᵢ + Ca∞ )])

# from the reported area value in Prinz 2003; this is actually the area of a cylinder
# without the ends (2πrh) with dimensions r = 25µm, h = 400µm as reported by Liu 
# here we just get a fictive 
radius = sqrt(0.628e-3)/pi * cm

@named neuron = Soma([NaV, CaS, CaT, KA, KCa, Kdr, H, leak], gradients,
                     radius = radius, V0 = -50mV, applied = 0nA, aux = [calcium_conversion])

t = 3000
sim = Simulation(neuron, time = t*ms)

solution = solve(sim, Rosenbrock23())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(t*5))

