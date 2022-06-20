# Compartments

## CompartmentSystem

A single `CompartmentSystem` is sufficient to model a simple neuron--either as an
isopotential sphere- or cylindrical-shaped neuron, or as a "point" neuron with no morphological geometry. At
minimum, a `CompartmentSystem` will contain information about its instrinsic membrane
currents and ionic gradients (i.e. equilbrium potentials). We define these by directly
providing vectors of `ConductanceSystem` and `EquilibriumPotential` to the
`CompartmentSystem`:

```jldoctest compartment_doctest
using Conductor, Unitful, ModelingToolkit, IfElse
import Unitful: mV, mS, cm

Vₘ = MembranePotential()

nav_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

kdr_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4, name = :n)]

# Note: `IonChannel` is a convenience constructor that returns a `ConductanceSystem`
@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Sodium => 50.0mV, Potassium => -77.0mV, Leak => -54.4mV])

@named neuron = CompartmentSystem(Vₘ, channels, reversals)

# output
Model neuron with 10 equations
States (12):
  Ileak(t)
  Vₘ(t) [defaults to -60.0]
  IKdr(t)
  INaV(t)
  NaV₊Vₘ(t)
  NaV₊g(t)
  NaV₊h(t)
  NaV₊m(t)
⋮
Parameters (8):
  EK
  El
  cₘ
  aₘ
  ENa
  NaV₊gbar
  leak₊gbar
  Kdr₊gbar
```

!!! note
    We use `IfElse.ifelse` to handle discontinuities in the rate equations for Sodium and
    Potassium channel dynamics. In Julia v1.8+ `Base.ifelse` may be used directly without
    the need for an external dependency. 

In the case above, the `CompartmentSystem` constructor assumes a dimensionless `geometry =
Point()`. The maximum magnitudes of the ion channel conductances, ``\\overline{g}``, have
units of `SpecificConductance` (mS/cm²) that scale with the surface area of the compartment.
However, when the geometry is a `Point()`, the `CompartmentSystem` ignores surface area
scaling (internally, the area is fixed to 1.0). We can instead provide our own `<:Geometry`
object to describe the shape and size of the compartment:

```jldoctest compartment_doctest; output = false
import Unitful: µm
soma_shape = Sphere(radius = 20µm)

@named neuron = CompartmentSystem(Vₘ, channels, reversals; geometry = soma_shape)

# output
Model neuron with 10 equations
States (12):
  Ileak(t)
  Vₘ(t) [defaults to -60.0]
  IKdr(t)
  INaV(t)
  NaV₊Vₘ(t)
  NaV₊g(t)
  NaV₊h(t)
  NaV₊m(t)
⋮
Parameters (8):
  EK
  El
  cₘ
  aₘ
  ENa
  NaV₊gbar
  leak₊gbar
  Kdr₊gbar
```

If we try to simulate the neuron we've modeled so far, the result isn't very interesting:

```@setup compartment_example
using Conductor, Unitful, ModelingToolkit, IfElse
import Unitful: mV, mS, cm, µm
Vₘ = MembranePotential()
nav_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]
kdr_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4, name = :n)]
@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)
channels = [NaV, Kdr, leak];
reversals = Equilibria([Sodium => 50.0mV, Potassium => -77.0mV, Leak => -54.4mV])
soma_shape = Sphere(radius = 20µm)
@named neuron = CompartmentSystem(Vₘ, channels, reversals; geometry = soma_shape)
```

```@example compartment_example
import Unitful: ms
using OrdinaryDiffEq, Plots
sim = Simulation(neuron, time = 300ms)
solution = solve(sim, Rosenbrock23())
plot(solution; vars=[Vₘ])
savefig("silentneuron_plot.svg"); nothing # hide
```
![](silentneuron_plot.svg)

The neuron isn't spontaneously active. To make the neuron produce spikes, we can write an
equation for an electrode current and provide it to `CompartmentSystem`: 

```jldoctest compartment_doctest
import Unitful: µA, pA
@named Iₑ = IonCurrent(NonIonic)

# A 400 picoamp squarewave pulse when 100ms > t > 200ms
electrode_pulse = Iₑ ~ IfElse.ifelse(t > 100.0,
                                     IfElse.ifelse(t < 200.0,
                                                   ustrip(Float64, µA, 400pA),
                                                   0.0),
                                     0.0)

@named neuron_stim = CompartmentSystem(Vₘ, channels, reversals;
                                       geometry = soma_shape,
                                       stimuli = [electrode_pulse])
# output
Model neuron_stim with 11 equations
States (13):
  Ileak(t)
  Iₑ(t) [defaults to ifelse(t > 100.0, ifelse(t < 200.0, 0.0004, 0.0), 0.0)]
  Vₘ(t) [defaults to -60.0]
  IKdr(t)
  INaV(t)
  NaV₊Vₘ(t)
  NaV₊g(t)
  NaV₊h(t)
⋮
Parameters (8):
  EK
  El
  cₘ
  aₘ
  ENa
  NaV₊gbar
  leak₊gbar
  Kdr₊gbar
```

Putting it all together, our neuron simulation now produces a train of action potentials:

```@setup compartment_example
import Unitful: µA, pA
@named Iₑ = IonCurrent(NonIonic)

# A 400 picoamp squarewave pulse when 100ms > t > 200ms
electrode_pulse = Iₑ ~ IfElse.ifelse(t > 100.0,
                                     IfElse.ifelse(t < 200.0,
                                                   ustrip(Float64, µA, 400pA),
                                                   0.0),
                                     0.0)

@named neuron_stim = CompartmentSystem(Vₘ, channels, reversals;
                                  geometry = soma_shape,
                                  stimuli = [electrode_pulse])
```

```@example compartment_example
sim = Simulation(neuron_stim, time = 300ms)
solution = solve(sim, Rosenbrock23())
plot(solution; vars=[Vₘ])
savefig("spikingneuron_plot.svg"); nothing # hide
```
![](spikingneuron_plot.svg)

## MultiCompartmentSystem

```@docs
CompartmentSystem
Junction
MultiCompartmentSystem
MultiCompartment
```
