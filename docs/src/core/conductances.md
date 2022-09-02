# [ConductanceSystem](@id conductance)

A `ConductanceSystem` describes the state- and time-dependence of a conductance (i.e.
classically _g_, the inverse of resistance; typically measured in Siemens). It can be used
to model the conductance associated with ionic membrane currents, synaptic currents, and
axial currents that flow between connected neuronal compartments. By default, a `Conductance
System` is composed of zero or more gates and a parameter for the maximum conductance value
(``\overline{g}``). The output of a `ConductanceSystem` is equal to the the product of each
individual gate and ``\overline{g}``.

To model sodium conductance in a Hodgkin-Huxley model:

```@example
using Conductor, ModelingToolkit, Unitful
import Unitful: mV, mS, cm

Vₘ = MembranePotential()

nav_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@assert length(equations(NaV)) == 3; # hide

equations(NaV) # includes: g(t) ~ gbar*(m(t)^3)*h(t)
```

```@docs
ConductanceSystem
IonChannel
AxialConductance
SynapticChannel
```
