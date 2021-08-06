cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__); Pkg.instantiate()

include("hodgkin_huxley_test.jl")

byhand_prob = ODEProblem{true}(hodgkin_huxley!, u0, (0.,300.), p)
byhand_sol = solve(byhand_prob, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);
current_mtk_sol = solve(sim, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);

tsteps = 0.0:0.025:300.0
byhand_out = Array(byhand_sol(tsteps, idxs=1))
current_mtk_out = Array(current_mtk_sol(tsteps, idxs=2))

@test isapprox(byhand_out, current_mtk_out, rtol=0.0001)
