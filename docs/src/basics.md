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

Each conductance is associated with a [`CompartmentSystem`](@ref compartment). The
properties of the parent compartment (for example, equilibrium potentials, membrane capacitance,
and geometry), influence the magnitude of intracellular current produced by conductances.

Finally, connections between neurons can be explicitly defined with [`Synapse`](@ref
synapses), which
defines a directed edge from a presynaptic compartment to a postsynaptic compartment. A
[`NeuronalNetworkSystem`](@ref networks) constructs and manages the equations that govern synaptic
conductances between neurons.

## Simple Two Neuron Simulation

Here's a quick copy-paste example of two synaptically coupled neurons with Hodgkin-Huxley
dynamics:

```julia
using Conductor, IfElse, OrdinaryDiffEq, Plots, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K

Vₘ = MembranePotential(-65mV)

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
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named Iₑ = IonCurrent(NonIonic)
holding_current = Iₑ ~ ustrip(Float64, µA, 5000pA)

@named neuron1 = Compartment(Vₘ, channels, reversals;
                             geometry = Cylinder(radius = 25µm, height = 400µm),
                             stimuli = [holding_current])

@named neuron2 = Compartment(Vₘ, channels, reversals;
                             geometry = Cylinder(radius = 25µm, height = 400µm))
                                   
# Synaptic model
Vₓ = ExtrinsicPotential()
syn∞ = 1/(1 + exp((-35 - Vₓ)/5))
τsyn = (1 - syn∞)/(1/40)
syn_kinetics = Gate(SteadyStateTau, syn∞, τsyn, name = :z)
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@named Glut = SynapticChannel(Cation, [syn_kinetics]; max_s = 30nS);

net = NeuronalNetworkSystem([Synapse(neuron1 => neuron2, Glut, EGlut)])

total_time = 250
sim = Simulation(net, time = total_time*ms)

solution = solve(sim, Rosenbrock23())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(total_time*5))

```

