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



