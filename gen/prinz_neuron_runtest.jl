cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__); Pkg.instantiate()

using OrdinaryDiffEq, Test

include("prinz_neuron_setup.jl")

byhand_prob = ODEProblem{true}(prinz_neuron!, u0, (0.,2000.), p)
byhand_sol = solve(byhand_prob, Rosenbrock23(), reltol=1e-8, abstol=1e-8, saveat=0.025);
current_mtk_sol = solve(sim, Rosenbrock23(), reltol=1e-8, abstol=1e-8, saveat=0.025);

tsteps = 0.0:0.025:2000.0
byhand_out = Array(byhand_sol(tsteps, idxs=2))
current_mtk_out = current_mtk_sol[neuron.sys.Vₘ]

@test isapprox(byhand_out, current_mtk_out, rtol=0.001)
