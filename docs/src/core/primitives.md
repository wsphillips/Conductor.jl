# Constants & Primitives

## [Primitives](@id primitives)

Primitives are special symbolic variables that are annotated with additional domain-specific
metadata. Conductor.jl uses primitives found in user-defined equations to validate models
and improve code generation.

```@docs
MembranePotential
ExtrinsicPotential
IonConcentration
IonCurrent
EquilibriumPotential
Conductor.Temperature
```

## Constants

```@docs
Conductor.t
Conductor.D
Conductor.ℱ
Conductor.IonSpecies
Conductor.R
Conductor.mR
```
### Quantities

In addition, Conductor.jl defines the following useful Unitful.jl-compatible units:

```@docs
Conductor.SpecificConductance
Conductor.SpecificCapacitance
```
