# Conductor.jl

**Conductor.jl is a WIP. If the idea of a Julia-based neuronal network simulator engine sounds exciting to you, please feel free to reach out** 

![GitHub Workflow Status (branch)](https://img.shields.io/github/workflow/status/wsphillips/Conductor.jl/CI?label=tests&logo=julia)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://wsphillips.github.io/Conductor.jl/stable)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://wsphillips.github.io/Conductor.jl/dev)

Conductor.jl aims to be a platform for quickly and flexibly building high-performance,
multi-scale neuronal network models in Julia. Under the hood it's being built on top of
ModelingToolkit.jl--so all the tools available in the SciML and DiffEq ecosystem are (or
soon will be) useable and composable with the neuronal models built here.

To install, tagged releases are available through the public registry:

```julia
# From Julia REPL
]add Conductor
```

While Conductor.jl is still in early development, you can get a feel for what's going on by looking in
the `demo` directory of this repository. Clone the latest tagged version of this repository:

```
git clone https://github.com/wsphillips/Conductor.jl
cd Conductor.jl
git checkout v0.0.2
```
Then from a Julia REPL:
```julia
cd("/path/to/Conductor.jl/demo")
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

You should then be able to open and step through the various demo script examples.

## Acknowldegements

Conductor.jl is based on the acausal component modeling paradigm in
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). The initial draft of
Conductor.jl was derived from an implementation of a stomatogastric ganglion (STG) model,
which was written in Julia by [Dhruva Raman](http://www-control.eng.cam.ac.uk/Main/DhruvaRaman), and based on published works by [Astrid Prinz](http://www.biology.emory.edu/research/Prinz/) et
al.

The original Julia/ModelingToolkit STG model template: [NeuronBuilder.jl](https://github.com/Dhruva2/NeuronBuilder)

STG model papers:

Prinz et al. 2003 *The functional consequences of changes in the strength and duration of synaptic inputs to
oscillatory neurons* [J. Neuroscience](https://www.jneurosci.org/content/23/3/943.full)

Prinz et al. 2003 *Alternative to hand-tuning conductance-based models: construction and
analysis of databases of model neurons* [J.
Neurophysiology](https://journals.physiology.org/doi/full/10.1152/jn.00641.2003)

Prinz et al. 2004 *Similar network activity from disparate circuit parameters* [Nature
Neuroscience](https://www.nature.com/articles/nn1352)

Thanks also to [Srinivas Gorur-Shandilya](https://srinivas.gs/) for advice and
contributions related to model implementation.
