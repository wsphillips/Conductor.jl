# Basics

## Overview

Conductor.jl builds neurobiological models in a bottom-up workflow, beginning with
user-defined equations. Equations may be written using [primitives](@ref primitives)--declared
symbolic variables that represent physical quantities like voltage, current, ion
concentrations, etc.

State dynamics in neurons are described predominately by [gates](@ref gates), which convert
arbitrary inputs into activation weights. A [`ConductanceSystem`](@ref conductance) is
comprised of one or more gates. Conductances combine the outputs of gates, producing a
scalar quantity that can represent the permeability of ion channels, synaptic channels, and
axial conductances between morphologically connected compartments.

Each conductance is associated with a [`CompartmentSystem`](@ref compartments). The
properties of the parent compartment (for example, equilibrium potentials, membrane capacitance,
and geometry) influence the magnitude of intracellular current produced by conductances.

Finally, connections between neurons can be explicitly defined with a `NetworkTopology`,
which stores a multilayer graph representation of the neuronal network. A
[`NeuronalNetworkSystem`](@ref networks) constructs and manages the equations that govern
synaptic conductances between neurons.

## Simple Two Neuron Simulation

Here's a quick copy-paste example of two synaptically coupled neurons with Hodgkin-Huxley
dynamics:

```@example
using Conductor,  OrdinaryDiffEq, Plots, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K

Vₘ = MembranePotential(-65mV)

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

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named Iₑ = IonCurrent(NonIonic)
@named I_hold = IonCurrent(NonIonic, 5000pA, dynamic = false)
holding_current = Iₑ ~ I_hold

dynamics_1 = HodgkinHuxley(Vₘ, channels, reversals;
                         geometry = Cylinder(radius = 25µm, height = 400µm),
                         stimuli = [holding_current])
dynamics_2 = HodgkinHuxley(Vₘ, channels, reversals;
                           geometry = Cylinder(radius = 25µm, height = 400µm))

@named neuron1 = Compartment(dynamics_1)
@named neuron2 = Compartment(dynamics_2)

# Synaptic model
Vₓ = ExtrinsicPotential()
syn∞ = 1/(1 + exp((-35 - Vₓ)/5))
τsyn = (1 - syn∞)/(1/40)
syn_kinetics = Gate(SteadyStateTau, syn∞, τsyn, name = :z)
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@named Glut = SynapticChannel(Cation, [syn_kinetics]; max_s = 30nS);

topology = NetworkTopology([neuron1, neuron2], [Glut]);
topology[neuron1, neuron2] = Glut
reversal_map = Dict([Glut => EGlut])

@named net = NeuronalNetworkSystem(topology, reversal_map)
total_time = 250
sim = Simulation(net, time = total_time*ms)

solution = solve(sim, Rosenbrock23(), abstol=1e-3, reltol=1e-3, saveat=0.2)
plot(solution; vars=[neuron1.Vₘ, neuron2.Vₘ])
```

## Step-by-step explanation
```@setup gate_example
using Conductor,  OrdinaryDiffEq, Plots, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K
```
We start by defining a primitive variable for voltage, `Vₘ`, which we'll use to write
equations describing the kinetics of some Hodgkin-Huxley gating particles.

```@example gate_example; continued = true
Vₘ = MembranePotential(-65mV)
```
Often the voltage- and time-dependent dynamics of ion channels are described in books and
journal articles using forward, ``\alpha``, and reverse, ``\beta``, reaction rates:

```math
\begin{aligned}
\alpha(V_{m}) &= 0.1(V_{m}+40)/1-e^{-(V_{m}+40)/10} \\
\beta(V_{m}) &= 4e^{-(V_{m}+65)/18}
\end{aligned}
```
To reproduce the math, we'll define a `Gate` and indicate the format of our equations with a
trait, `AlphaBeta`. We'll also use an optional keyword argument, `p`, which sets the
exponent of the resulting gating variable (for example, `p = 3` to set ``m^{3}`` as found in
Hodgkin-Huxley-style Sodium channels). 

```@example gate_example; continued = true
# Fast Sodium channel kinetics
nav_kinetics = [
    # the activation gate, m
    Gate(AlphaBeta,
         ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    # the inactivation gate, h
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

# Delayed rectifier Potassium kinetics
kdr_kinetics = [
    # the activation gate, n
    Gate(AlphaBeta,
         ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4, name = :n)]
```
!!! note
    There's a discontinuity in the original equations, so we use `ifelse` to avoid a
    divide-by-zero error.

Now that we've defined the gating variables, we can construct some ion channels and define a
set of matching primitives representing equilibrium potentials for sodium, potassium, and
non-specific leak current. 

```@example gate_example; continued = true
@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)
channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])
```
We also need to model current injection to stimulate spiking. We'll declare a new primitive
current and use the resulting symbolic to write an equation describing what the value of the
electrode current should be.

```@example gate_example; continued = true
@named Iₑ = IonCurrent(NonIonic)
@named I_hold = IonCurrent(NonIonic, 5000pA, dynamic = false)
holding_current = Iₑ ~ I_hold
```

Finally, we construct our two neurons, providing the holding current stimulus to `neuron1`,
which will be our presynaptic neuron.

```@example gate_example
dynamics_1 = HodgkinHuxley(Vₘ, channels, reversals;
                         geometry = Cylinder(radius = 25µm, height = 400µm),
                         stimuli = [holding_current])

dynamics_2 = HodgkinHuxley(Vₘ, channels, reversals;
                           geometry = Cylinder(radius = 25µm, height = 400µm))

@named neuron1 = Compartment(dynamics_1)
@named neuron2 = Compartment(dynamics_2)
``` 
For our neurons to talk to each other, we'll need a model for a synaptic conductance. This
time we'll use a model that's presented in a different form in the literature. A model of a
glutamatergic excitatory synapse adapted from
[Prinz et al 2004](https://www.nature.com/articles/nn1352):

```math
    \begin{aligned}
    s_{\infty}(V_{pre} &= \frac{1}{1 + e^{-35 - V_{pre}/5}} \\
    \tau_{s}   &= \frac{1 - s_{\infty}(V_{pre}}{1/40}
    \end{aligned}
```
We will start by declaring a primitive, `ExtrinsicPotential`, which represents a voltage
coming from an outside source (in this case, a presynaptic neuron). 

```@example gate_example; continued=true
Vₓ = ExtrinsicPotential()
```
Next, we define a `Gate` and indicate that our kinetic equations have the form
`SteadyStateTau`. 

```@example gate_example; continued=true
syn∞ = 1/(1 + exp((-35 - Vₓ)/5))
τsyn = (1 - syn∞)/(1/40)
syn_kinetics = Gate(SteadyStateTau, syn∞, τsyn, name = :z)
```
We'll then build a `SynapticChannel` and define an equilibrium potential.

```@example gate_example; continued=true
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@named Glut = SynapticChannel(Cation, [syn_kinetics]; max_s = 30nS);
```
Given our two neurons and a synaptic channel, we can model a miniature circuit by defining
a synapse between the neurons and then constructing a `NeuronalNetworkSystem`:

```@example gate_example; continued=true
topology = NetworkTopology([neuron1, neuron2], [Glut]);
add_synapse!(topology, neuron1, neuron2, Glut)
reversal_map = Dict([Glut => EGlut])

@named net = NeuronalNetworkSystem(topology, reversal_map)
```
Now we're ready to run our simulation.

```@example gate_example
total_time = 250
sim = Simulation(net, time = total_time*ms)
solution = solve(sim, Rosenbrock23(), abstol=1e-3, reltol=1e-3, saveat=0.2)
plot(solution; vars=[neuron1.Vₘ, neuron2.Vₘ])
```
