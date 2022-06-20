# Compartments

## CompartmentSystem

A single `CompartmentSystem` is sufficient to model a simple neuron--either as an
isopotential sphere- or cylindrical-shaped neuron, or as a "point" neuron with no morphological geometry. At
minimum, a `CompartmentSystem` will contain information about its instrinsic membrane
currents and ionic gradients (i.e. equilbrium potentials). We define these by directly
providing vectors of `ConductanceSystem` and `EquilibriumPotential` to the
`CompartmentSystem`:

```@example compartment_example
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
@assert length.((equations(neuron), states(neuron), parameters(neuron))) == (10,12,8) # hide
neuron # hide
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

```@example compartment_example
import Unitful: µm
soma_shape = Sphere(radius = 20µm)

@named neuron = CompartmentSystem(Vₘ, channels, reversals; geometry = soma_shape)
@assert length.((equations(neuron), states(neuron), parameters(neuron))) == (10,12,8); # hide
nothing # hide
```

If we try to simulate the neuron we've modeled so far, the result isn't very interesting:

```@example compartment_example
import Unitful: ms
using OrdinaryDiffEq, Plots
sim = Simulation(neuron, time = 300ms)
solution = solve(sim, Rosenbrock23())
plot(solution; vars=[Vₘ])
```
The neuron isn't spontaneously active. To make the neuron produce spikes, we can write an
equation for an electrode current and provide it to `CompartmentSystem`: 

```@example compartment_example
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
@assert length.((equations(neuron_stim), states(neuron_stim), parameters(neuron_stim))) == (11,13,8); # hide
neuron_stim # hide
```
Putting it all together, our neuron simulation now produces a train of action potentials:

```@example compartment_example
sim = Simulation(neuron_stim, time = 300ms)
solution = solve(sim, Rosenbrock23())
plot(solution; vars=[Vₘ])
```

## MultiCompartmentSystem

Neurons with multiple compartments are explicitly constructed through the use of `Junction`,
which denotes an connection (edge) between two compartments (nodes). By default, a
`Junction` is assumed to be symmetrical--the axial conductance from
compartment A to B is the same as from compartment B to A. If conductance between
compartments is asymmetric, two `Junction` objects must be defined--one for conductance in
each direction.





```@docs
CompartmentSystem
Junction
MultiCompartmentSystem
MultiCompartment
```
