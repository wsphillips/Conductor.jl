using Conductor, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, ms, mM, µM, µA
import Conductor: Na, K, Ca, Cation, Leak

Vₘ = MembranePotential(-50mV)
Caᵢ = Concentration(Calcium, 0.05µM, dynamic = true)
ICa = IonCurrent(Calcium, aggregate = true)

nav_kinetics = [
    Gate(SteadyStateTau,
         1.0/(1.0 + exp((Vₘ + 25.5)/-5.29)),
         2.64 - (2.52/(1 + exp((Vₘ + 120.0)/-25.0))),
         p = 3, name = :m)
    Gate(SteadyStateTau,
         1.0/(1.0 + exp((Vₘ + 48.9)/5.18)),
         (1.34/(1.0 + exp((Vₘ + 62.9)/-10.0)))*(1.5 + 1.0/(1.0 + exp((Vₘ + 34.9)/3.6))), name = :h)]
         
cas_kinetics = [
    Gate(SteadyStateTau,
         1.0/(1.0+exp((Vₘ+33.0)/-8.1)),
         2.8 + 14.0/(exp((Vₘ+27.0)/10.0) + exp((Vₘ+70.0)/-13.0)),
         p = 3, name = :m)
    Gate(SteadyStateTau,
         1.0/(1.0+exp((Vₘ+60.0)/6.2)),
         120.0 + 300.0/(exp((Vₘ+55.0)/9.0) + exp((Vₘ+65.0)/-16.0)), name = :h)]

cat_kinetics = [
    Gate(SteadyStateTau,
         1.0/(1.0 + exp((Vₘ+27.1)/-7.2)),
         43.4 - 42.6/(1.0 + exp((Vₘ+68.1)/-20.5)),
         p = 3, name = :m)
    Gate(SteadyStateTau,
         1.0/(1.0 + exp((Vₘ+32.1)/5.5)),
         210.0 - 179.6/(1.0 + exp((Vₘ+55.0)/-16.9)), name = :h)]

ka_kinetics = [
    Gate(SteadyStateTau,
         1.0/(1.0+exp((Vₘ+27.2)/-8.7)),
         23.2 - 20.8/(1.0+exp((Vₘ+32.9)/-15.2)),
         p = 3, name = :m)
    Gate(SteadyStateTau,
         1.0/(1.0+exp((Vₘ+56.9)/4.9)),
         77.2 - 58.4/(1.0+exp((Vₘ+38.9)/-26.5)), name = :h)]

kca_kinetics = [
    Gate(SteadyStateTau,
         (Caᵢ/(Caᵢ+3.0))/(1.0+exp((Vₘ+28.3)/-12.6)),
         180.6 - 150.2/(1.0+exp((Vₘ+46.0)/-22.7)),
         p = 4, name = :m)]
       
kdr_kinetics = [
    Gate(SteadyStateTau,
         1.0/(1.0+exp((Vₘ+12.3)/-11.8)),
         14.4 - 12.8/(1.0+exp((Vₘ+28.3)/-19.2)),
         p = 4, name = :m)]

h_kinetics = [
    Gate(SteadyStateTau,
         1.0/(1.0+exp((Vₘ+75.0)/5.5)),
         2/( exp((Vₘ+169.7)/(-11.6)) + exp((Vₘ- 26.7)/(14.3))), name = :m)]

Ca_Nernst(Ca) = 1000*((-8.314*283.15)/(2 * 96485.365))*log(max(Ca,0.001)/3000.0)
@register_symbolic Ca_Nernst(Ca)
ModelingToolkit.get_unit(op::typeof(Ca_Nernst), args) = mV

gradients = Equilibria(Pair[Sodium    =>  50.0mV,
                            Potassium => -80.0mV,
                            Cation    => (-20mV, :H),
                            Leak      => -50mV,
                            # x1000 for mV; thresholding Cai in case negative conc
                            # TODO: write nernst function + use Unitful constants (RT/zℱ)
                            Calcium   => Ca_Nernst(Caᵢ)]);

# Reported area value in Prinz 2003; the area of a cylinder without the ends (2πrh).
# Dimensions r = 25µm, h = 400µm given by Liu et al 1998.

#r = 25µm; h = 400µm; area = round(ustrip(Float64, cm^2, 2π*r*h), sigdigits=3)

geo = Cylinder(radius = 25µm, height = 400µm)

# In the Prinz model, the f value of "14.96" appears to come from canceling out Cm in the
# Liu f value: 0.94/0.628nF * 10 = 14.96. There's no explanation for why it's bumped up by a
# factor of 10. The time constant being longer was explained in Goldman 2001 J. Neuroscience.
# Internal current units are µA (mS * mV) -> x1000 to match Prinz/Liu factor (uses nA)

# conversion as published
@parameters τCa = 200.0 [unit=ms] Ca∞ = 0.05 [unit=µM] fₚ = 14.96 [unit=µM/µA]
@named calcium_conversion = ODESystem([D(Caᵢ) ~ ((-fₚ*ICa*1000) - Caᵢ + Ca∞)/τCa]);

@named NaV = IonChannel(Sodium, nav_kinetics) 
@named CaT = IonChannel(Calcium, cat_kinetics)
@named CaS = IonChannel(Calcium, cas_kinetics)
@named KA  = IonChannel(Potassium, ka_kinetics)
@named KCa = IonChannel(Potassium, kca_kinetics)
@named Kdr = IonChannel(Potassium, kdr_kinetics)
@named H   = IonChannel(Cation, h_kinetics)
@named leak = IonChannel(Leak)

