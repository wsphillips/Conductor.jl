cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__); Pkg.instantiate()

include("synapse_setup.jl")

indexof(sym,syms) = findfirst(isequal(sym), syms)

simplified_system = structural_simplify(network)
# since the index of states is non-deterministic, we look it up at runtime
idx = indexof(neuron1.sys.Vₘ, states(simplified_system))

byhand_prob = ODEProblem{true}(simple_synapse!, u0, (0.,t), p)
byhand_sol = solve(byhand_prob, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);
current_mtk_sol = solve(simul, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);

tsteps = 0.0:0.025:t
byhand_out = Array(byhand_sol(tsteps, idxs=3))
current_mtk_out = Array(current_mtk_sol(tsteps, idxs=idx))

using Test

@test isapprox(byhand_out, current_mtk_out, rtol=0.0001)
