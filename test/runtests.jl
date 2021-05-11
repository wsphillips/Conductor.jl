
using Conductor
import Conductor: Gate, SteadyStateTau, IonChannel, Sodium
Vₘ = Conductor.MembranePotential

nav_kinetics = [
    Gate(SteadyStateTau();
         m∞ = 1.0/(1.0 + exp((Vₘ + 25.5)/-5.29)),
         τₘ = 1.32 - 1.26/(1 + exp((Vₘ + 120.0)/-25.0)), p = 3)
    Gate(SteadyStateTau();
         h∞ = 1.0/(1.0 + exp((Vₘ + 48.9)/5.18)),
         τₕ = (0.67/(1.0 + exp((Vₘ + 62.9)/-10.0)))*(1.5 + 1.0/(1.0 + exp((Vₘ + 34.9)/3.6))))]

@named NaV = IonChannel(Sodium, nav_kinetics, 120mS/cm^2) 

cas_kinetics = [
    Gate(m∞ = 1.0/(1.0+exp((V+33.0)/-8.1)),
         τₘ = 1.4 + 7.0/(exp((V+27.0)/10.0) + exp((V+70.0)/-13.0)), p = 3)
    Gate(h∞ = 1.0/(1.0 + exp((Vₘ + 48.9)/5.18)),
         τₕ = (0.67/(1.0 + exp((Vₘ + 62.9)/-10.0)))*(1.5 + 1.0/(1.0 + exp((Vₘ + 34.9)/3.6))))]

@named CaS = IonChannel(Calcium, cas_kinetics, 20mS/cm^2)
