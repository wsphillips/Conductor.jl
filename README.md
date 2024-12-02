# Conductor.jl

<p align="center">
<img src="https://img.shields.io/github/actions/workflow/status/wsphillips/Conductor.jl/CI.yml?branch=main&label=Tests&logo=julia&style=for-the-badge" alt="GitHub Workflow Status" />
<a href="https://wsphillips.github.io/Conductor.jl/stable">
<img src="https://img.shields.io/badge/docs-stable-blue.svg?style=for-the-badge" />
</a>
<a href="https://wsphillips.github.io/Conductor.jl/dev">
<img src="https://img.shields.io/badge/docs-dev-blue.svg?style=for-the-badge"/>
</p>

![LBFGS Optimization](lbfgs_stgneuron.gif)

Conductor.jl aims to be a platform for quickly and flexibly building high-performance,
multi-scale neuronal network models in Julia. Under the hood it's being built on top of
ModelingToolkit.jl--so all the tools available in the SciML and DiffEq ecosystem are (or
soon will be) useable and composable with the neuronal models built here.

To install, tagged releases are available through the public registry:

```julia
# From Julia REPL
]add Conductor
```

You can get a feel for what's going on by looking in
the `demo` directory of this repository. Clone the repository:

```
git clone https://github.com/wsphillips/Conductor.jl
```
Then from a Julia REPL:
```julia
cd("/path/to/Conductor.jl/demo")
using Pkg; Pkg.activate("."); Pkg.instantiate()
```

You should then be able to open and step through the various demo script examples.

## Acknowledgements

Conductor.jl is based on the acausal component modeling paradigm in
[ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl). The project began with an implementation of a stomatogastric ganglion (STG) model,
which was written in Julia by [Dhruva Raman](http://www-control.eng.cam.ac.uk/Main/DhruvaRaman), and based on published works by [Astrid Prinz](http://www.biology.emory.edu/research/Prinz/) et
al.

The original Julia STG model: [NeuronBuilder.jl](https://github.com/Dhruva2/NeuronBuilder)

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
