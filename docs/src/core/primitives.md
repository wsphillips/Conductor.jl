# Constants & Primitives

# Primitives

Primitives are special symbolic variables that are annotated with additional domain-specific
metadata. Conductor.jl uses primitives found in user-defined equations to validate models
and improve code generation.

```@docs
MembranePotential
ExtrinsicPotential
IonConcentration
IonCurrent
EquilibriumPotential
```

# Constants

*IonSpecies* - an enumeration of common ions
    * Sodium
    * Potassium
    * Calcium
    * Chloride
