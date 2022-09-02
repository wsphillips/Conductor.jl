using Conductor, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, ms, mM, µM, µA
import Conductor: Na, K, Ca, Cation, Leak

Vₘ = MembranePotential(-50mV)
Caᵢ = Concentration(Calcium, 0.05µM, dynamic = true)
ICa = IonCurrent(Calcium, aggregate = true)

NaV₊m∞(V)    = 1.0/(1.0 + exp((V + 25.5)/-5.29))
NaV₊τₘ(V)    = 2.64 - (2.52/(1 + exp((V + 120.0)/-25.0)))
NaV₊h∞(V)    = 1.0/(1.0 + exp((V + 48.9)/5.18))
NaV₊τₕ(V)    = (1.34/(1.0 + exp((V + 62.9)/-10.0)))*(1.5 + 1.0/(1.0 + exp((V + 34.9)/3.6)))

CaS₊m∞(V)    = 1.0/(1.0+exp((V+33.0)/-8.1))
CaS₊τₘ(V)    = 2.8 + 14.0/(exp((V+27.0)/10.0) + exp((V+70.0)/-13.0))
CaS₊h∞(V)    = 1.0/(1.0+exp((V+60.0)/6.2))
CaS₊τₕ(V)    = 120.0 + 300.0/(exp((V+55.0)/9.0) + exp((V+65.0)/-16.0))

CaT₊m∞(V)    = 1.0/(1.0 + exp((V+27.1)/-7.2))
CaT₊τₘ(V)    = 43.4 - 42.6/(1.0 + exp((V+68.1)/-20.5))
CaT₊h∞(V)    = 1.0/(1.0 + exp((V+32.1)/5.5))
CaT₊τₕ(V)    = 210.0 - 179.6/(1.0 + exp((V+55.0)/-16.9))

KA₊m∞(V)     = 1.0/(1.0+exp((V+27.2)/-8.7))
KA₊τₘ(V)     = 23.2 - 20.8/(1.0+exp((V+32.9)/-15.2))
KA₊h∞(V)     = 1.0/(1.0+exp((V+56.9)/4.9))
KA₊τₕ(V)     = 77.2 - 58.4/(1.0+exp((V+38.9)/-26.5))

KCa₊m∞(V,Ca) = (Ca/(Ca + 3.0))/(1.0 + exp((V + 28.3)/-12.6))
KCa₊τₘ(V,Ca) = 180.6 - 150.2/(1.0 + exp((V + 46.0)/-22.7)) 

Kdr₊m∞(V)    = 1.0/(1.0+exp((V+12.3)/-11.8))
Kdr₊τₘ(V)    = 14.4 - 12.8/(1.0+exp((V+28.3)/-19.2))

H₊m∞(V)      = 1.0/(1.0+exp((V+75.0)/5.5))
H₊τₘ(V)      = 2/( exp((V+169.7)/(-11.6)) + exp((V- 26.7)/(14.3)))

@register_symbolic NaV₊m∞(V)   
@register_symbolic NaV₊τₘ(V)   
@register_symbolic NaV₊h∞(V)   
@register_symbolic NaV₊τₕ(V)   
@register_symbolic CaS₊m∞(V)   
@register_symbolic CaS₊τₘ(V)   
@register_symbolic CaS₊h∞(V)   
@register_symbolic CaS₊τₕ(V)   
@register_symbolic CaT₊m∞(V)   
@register_symbolic CaT₊τₘ(V)   
@register_symbolic CaT₊h∞(V)   
@register_symbolic CaT₊τₕ(V)   
@register_symbolic KA₊m∞(V)    
@register_symbolic KA₊τₘ(V)    
@register_symbolic KA₊h∞(V)    
@register_symbolic KA₊τₕ(V)    
@register_symbolic KCa₊m∞(V,Ca)
@register_symbolic KCa₊τₘ(V,Ca)
@register_symbolic Kdr₊m∞(V)   
@register_symbolic Kdr₊τₘ(V)   
@register_symbolic H₊m∞(V)     
@register_symbolic H₊τₘ(V)     

ModelingToolkit.get_unit(op::typeof(NaV₊m∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(NaV₊τₘ), args) = ms
ModelingToolkit.get_unit(op::typeof(NaV₊h∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(NaV₊τₕ), args) = ms
ModelingToolkit.get_unit(op::typeof(CaS₊m∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(CaS₊τₘ), args) = ms
ModelingToolkit.get_unit(op::typeof(CaS₊h∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(CaS₊τₕ), args) = ms
ModelingToolkit.get_unit(op::typeof(CaT₊m∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(CaT₊τₘ), args) = ms
ModelingToolkit.get_unit(op::typeof(CaT₊h∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(CaT₊τₕ), args) = ms
ModelingToolkit.get_unit(op::typeof(KA₊m∞) , args) = NoUnits
ModelingToolkit.get_unit(op::typeof(KA₊τₘ) , args) = ms
ModelingToolkit.get_unit(op::typeof(KA₊h∞) , args) = NoUnits
ModelingToolkit.get_unit(op::typeof(KA₊τₕ) , args) = ms
ModelingToolkit.get_unit(op::typeof(KCa₊m∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(KCa₊τₘ), args) = ms
ModelingToolkit.get_unit(op::typeof(Kdr₊m∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(Kdr₊τₘ), args) = ms
ModelingToolkit.get_unit(op::typeof(H₊m∞)  , args) = NoUnits
ModelingToolkit.get_unit(op::typeof(H₊τₘ)  , args) = ms

nav_kinetics = [
    Gate(SteadyStateTau, NaV₊m∞(Vₘ), NaV₊τₘ(Vₘ), p = 3, name = :m)
    Gate(SteadyStateTau, NaV₊h∞(Vₘ), NaV₊τₕ(Vₘ), name = :h)]
cas_kinetics = [
    Gate(SteadyStateTau, CaS₊m∞(Vₘ), CaS₊τₘ(Vₘ), p = 3, name = :m)
    Gate(SteadyStateTau, CaS₊h∞(Vₘ), CaS₊τₕ(Vₘ), name = :h)]

cat_kinetics = [
    Gate(SteadyStateTau, CaT₊m∞(Vₘ), CaT₊τₘ(Vₘ), p = 3, name = :m)
    Gate(SteadyStateTau, CaT₊h∞(Vₘ), CaT₊τₕ(Vₘ), name = :h)]

ka_kinetics = [
    Gate(SteadyStateTau, KA₊m∞(Vₘ), KA₊τₘ(Vₘ), p = 3, name = :m)
    Gate(SteadyStateTau, KA₊h∞(Vₘ), KA₊τₕ(Vₘ), name = :h)]

kca_kinetics = [
    Gate(SteadyStateTau, KCa₊m∞(Vₘ, Caᵢ), KCa₊τₘ(Vₘ, Caᵢ), p = 4, name = :m)]
       
kdr_kinetics = [
    Gate(SteadyStateTau, Kdr₊m∞(Vₘ), Kdr₊τₘ(Vₘ), p = 4, name = :m)]

h_kinetics = [ 
    Gate(SteadyStateTau, H₊m∞(Vₘ), H₊τₘ(Vₘ), name = :m)]

Ca_Nernst(Ca) = 1000*((-8.314*283.15)/(2 * 96485.365))*log(max(Ca,0.001)/3000.0)
@register_symbolic Ca_Nernst(Ca)

ModelingToolkit.get_unit(op::typeof(Ca_Nernst), args) = mV


