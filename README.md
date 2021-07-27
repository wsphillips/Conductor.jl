# Conductor.jl

**Conductor.jl is a WIP. If the idea of a Julia-based neuronal network simulator engine sounds exciting to you, please feel free to reach out** 

Conductor.jl aims to be a platform for quickly and flexibly building high-performance,
multi-scale neuronal network models in Julia. Under the hood it's being built on top of
ModelingToolkit.jl--so all the tools available in the SciML and DiffEq ecosystem are (or
soon will be) useable and composable with the neuronal models built here.

To install, I recommend cloning this repository to an easily reachable location, and adding
the directory via the package manager:

`git clone https://github.com/wsphillips/Conductor.jl`

```julia
# From Julia REPL
using Pkg; Pkg.add(path="/path/to/Conductor.jl")
```

While Conductor.jl is still in early development, you can get a feel for what's going on by looking in
the `demo` directory of this repository:

```julia
# From the Julia REPL
cd("/path/to/Conductor.jl/demo")
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

You should then be able to open and step through the various demo script examples.
