# [Compartments](@id compartments)

## Single Compartments

A `CompartmentSystem` is sufficient to model a simple neuron--either as an isopotential
spherical neuron, or as a "point" neuron with no morphological geometry. At minimum, a
`CompartmentSystem` will contain information about its intrinsic membrane currents and
ionic gradients (i.e. equilibrium potentials). We define these by directly providing vectors
of `ConductanceSystem` and `EquilibriumPotential` to the `CompartmentSystem`:

```@example compartment_example
using Conductor, Unitful, ModelingToolkit
import Unitful: mV, mS, cm

Vₘ = ParentScope(MembranePotential())

nav_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

kdr_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4, name = :n)]

# Note: `IonChannel` is a convenience constructor that returns a `ConductanceSystem`
@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Sodium => 50.0mV, Potassium => -77.0mV, Leak => -54.4mV])
dynamics = HodgkinHuxley(channels, reversals)
@named neuron = CompartmentSystem(Vₘ, dynamics)

@assert length.((equations(neuron), states(neuron), parameters(neuron))) == (10,10,8) # hide
neuron # hide
```

In the case above, the `CompartmentSystem` constructor assumes a dimensionless `geometry =
Point()`. The maximum magnitudes of the ion channel conductances, ``\overline{g}``, have
units of `SpecificConductance` (mS/cm²) that scale with the surface area of the compartment.
However, when the geometry is a `Point()`, the `CompartmentSystem` ignores surface area
scaling (internally, the area is fixed to 1.0). We can instead provide a `<:Geometry` object
to describe the shape and size of the compartment:

```@example compartment_example
import Unitful: µm
soma_shape = Sphere(radius = 20µm)
geo_dynamics = HodgkinHuxley(channels, reversals)
@named neuron = CompartmentSystem(Vₘ, geo_dynamics; geometry = soma_shape)
@assert length.((equations(neuron), states(neuron), parameters(neuron))) == (10,10,8); # hide
nothing # hide
```

If we try to simulate the neuron we've modeled so far, the result isn't very interesting:

```@example compartment_example
import Unitful: ms
using OrdinaryDiffEq, Plots
sim = Simulation(neuron, 300ms)
solution = solve(sim, Rosenbrock23())
plot(solution; idxs=[Vₘ])
```
The neuron isn't spontaneously active. To make the neuron produce spikes, we can write an
equation for an electrode current and provide it to `CompartmentSystem`: 

```@example compartment_example
import Unitful: µA, pA

# A 400 picoamp squarewave pulse when 100ms > t > 200ms
@named step_pulse = PulseTrain(amplitude = 400.0pA, duration = 100ms, delay = 100ms)

stim_dynamics = HodgkinHuxley(channels, reversals)

@named neuron_stim = CompartmentSystem(Vₘ, stim_dynamics;
                                       geometry = soma_shape,
                                       stimuli = [step_pulse])
@assert length.((equations(neuron_stim), states(neuron_stim), parameters(neuron_stim))) == (11,11,8); # hide
nothing # hide
```
Putting it all together, our neuron simulation now produces a train of action potentials:

```@example compartment_example
sim = Simulation(neuron_stim, 300ms)
solution = solve(sim, Rosenbrock23())
plot(solution; idxs=[Vₘ])
```

## Multiple Compartments

Neurons with multiple compartments are explicitly constructed by defining a
`MultiCompartmentTopology`, a directed graph that is provided to `MultiCompartmentSystem`.
Individual subcompartments are connected via `add_junction!`. By default, a junction between
compartments is assumed to be symmetrical--the axial conductance from compartment A to B is
the same value as from compartment B to A. However, if the conductance between compartments
is asymmetric, two junctions must be defined--one for the conductance in each direction.

For example usage, see the [`pinskyrinzel.jl`](https://github.com/wsphillips/Conductor.jl/blob/main/demo/pinskyrinzel.jl) demo.

```@docs
CompartmentSystem
MultiCompartmentSystem
MultiCompartment
```
