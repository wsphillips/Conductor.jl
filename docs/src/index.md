# Conductor.jl

## Introduction

Conductor.jl aims to be a platform for quickly and flexibly building high-performance,
multi-scale neuronal network models in Julia. Under the hood it's being built on top of
ModelingToolkit.jl--so all the tools available in the SciML and DiffEq ecosystem are (or
soon will be) useable and composable with the neuronal models built here.

## Installation

To install Conductor.jl, simply add it from the Pkg REPL:

```julia-repl
(@v1.7) pkg> add Conductor
```

Alternatively, you can use the following lines of Julia code:

```julia-repl
julia> using Pkg; Pkg.add(name="Conductor")
```
