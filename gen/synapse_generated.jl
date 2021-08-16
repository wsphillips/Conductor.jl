
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__); Pkg.instantiate()

include("synapse_setup.jl")
include("helpers.jl")

simplified = structural_simplify(network)

sim_exp = ModelingToolkit.build_torn_function(simplified; expression=true)
cleaned = clean_expr(sim_exp, simplified, rgf=true, pretty=true, name=Symbol("simple_synapse!"))

clipboard(format_text(string(cleaned), YASStyle()))


