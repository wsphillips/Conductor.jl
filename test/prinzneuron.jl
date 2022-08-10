
module PrinzSingleCompartment

using Test, OrdinaryDiffEq

include(normpath(@__DIR__, "..", "demo", "prinz_kinetics.jl"))

@testset "Prinz Single-Compartment Neuron" begin

channels = [NaV(100mS/cm^2),
            CaT(0mS/cm^2),
            CaS(4mS/cm^2),
            KA(0mS/cm^2),
            KCa(15mS/cm^2),
            Kdr(50mS/cm^2),
            H(.02mS/cm^2),
            leak(.03mS/cm^2)];

dynamics = HodgkinHuxley(Vₘ, channels, gradients; geometry = geo);
@named neuron = CompartmentSystem(dynamics, extensions = [calcium_conversion]);

@test length.([equations(neuron),
               states(neuron),
               parameters(neuron)]) == [39,39,17]

time = 2000 
simul_sys = Simulation(neuron, time = time*ms, return_system = true)

@test length.([equations(simul_sys),
               states(simul_sys),
               parameters(simul_sys)]) == [13,13,17]

# Prinz STG neuron hand-written reference implementation

# Steady-state functions
# Input vectors
V0 = -50.0
u0 = zeros(Float64, 13)
p  = zeros(Float64,17)

# Initial state
u0[1]  = 0.05 # Caᵢ
u0[2]  = V0 # Vₘ
u0[3]  = (1.0 + exp(-4.82042 - (0.189036V0)))^-1 # NaV₊m
u0[4]  = (1.0 + exp(9.44015 + 0.19305V0))^-1 # NaV₊h
u0[5]  = (1.0 + exp(-3.76389 - (0.138889V0)))^-1 # CaT₊m
u0[6]  = (1.0 + exp(5.83636 + 0.181818V0))^-1 # CaT₊h
u0[7]  = (1.0 + exp(-4.07407 - (0.123457V0)))^-1 # CaS₊m
u0[8]  = (1.0 + exp(9.67742 + 0.16129V0))^-1 # CaS₊h
u0[9]  = (1.0 + exp(-3.12644 - (0.114943V0)))^-1 # KA₊m
u0[10] = (1.0 + exp(11.6122 + 0.204082V0))^-1 # KA₊h
u0[11] = 0.05*((3.0 + 0.05)^-1)*((1.0 + exp(-2.24603 - (0.0793651V0)))^-1) # KCa₊m
u0[12] = (1.0 + exp(-1.04237 - (0.0847458V0)))^-1 # Kdr₊m
u0[13] = (1.0 + exp(13.6364 + 0.181818V0))^-1 # H₊m

# Parameters
p[1]  = 1.0 # cₘ      
p[2]  = area(geo) # aₘ      
p[3]  = 200 # τCa     
p[4]  = 0.05 # Ca∞     
p[5]  = 14.96 # f       
p[6]  = 50.0 # ENa     
p[7]  = -80.0 # EK      
p[8]  = -20.0 # EH      
p[9]  = -50.0 # El      
p[10] = 100.0 # NaV₊gbar
p[11] = 0.0 # CaT₊gbar  
p[12] = 4.0 # CaS₊gbar  
p[13] = 0.0 # KA₊gbar   
p[14] = 15.0 # KCa₊gbar 
p[15] = 50.0 # Kdr₊gbar 
p[16] = 0.02 # H₊gbar   
p[17] = 0.03 # leak₊g   

function prinz_neuron!(du, u, p, t)
    @inbounds let Caᵢ         = u[1],
                  Vₘ          = u[2],
                  NaV₊m       = u[3],
                  NaV₊h       = u[4],
                  CaT₊m       = u[5],
                  CaT₊h       = u[6],
                  CaS₊m       = u[7],
                  CaS₊h       = u[8],
                  KA₊m        = u[9],
                  KA₊h        = u[10],
                  KCa₊m       = u[11],
                  Kdr₊m       = u[12],
                  H₊m         = u[13],

                  cₘ          = p[1],
                  aₘ          = p[2],
                  τCa         = p[3],
                  Ca∞         = p[4],
                  f           = p[5],
                  ENa         = p[6],
                  EK          = p[7],
                  EH          = p[8],
                  El          = p[9],
                  NaV₊gbar    = p[10],
                  CaT₊gbar    = p[11],
                  CaS₊gbar    = p[12],
                  KA₊gbar     = p[13],
                  KCa₊gbar    = p[14],
                  Kdr₊gbar    = p[15],
                  H₊gbar      = p[16],
                  leak₊g      = p[17],
                  ECa = (Vₘ + 12.199306599503455 * log(0.0003333333333333333 * max(Caᵢ, 0.001)))
        let
            du[1] =
                (
                    Ca∞ - Caᵢ - 1000 * f *
                    (CaS₊gbar * aₘ * CaS₊h * CaS₊m^3 *
                    ECa +
                    CaT₊gbar * aₘ * CaT₊h * CaT₊m^3 *
                    ECa)
                   ) * inv(τCa)

            du[2] =
                (
                    -aₘ * leak₊g * (Vₘ - El) +
                    -H₊gbar * aₘ * H₊m * (Vₘ - EH) +
                    -KCa₊gbar * aₘ * KCa₊m^4 * (Vₘ - EK) +
                    -Kdr₊gbar * aₘ * Kdr₊m^4 * (Vₘ - EK) +
                    -CaS₊gbar * aₘ * CaS₊h * CaS₊m^3 * ECa +
                    -CaT₊gbar * aₘ * CaT₊h * CaT₊m^3 * ECa +
                    -KA₊gbar * aₘ * KA₊h * KA₊m^3 * (Vₘ - EK) +
                    -NaV₊gbar * aₘ * NaV₊h * NaV₊m^3 * (Vₘ - ENa)
                ) * inv(aₘ) * inv(cₘ)

            du[3] =
                (-NaV₊m + inv(1.0 + exp(-4.8204158790170135 - 0.1890359168241966Vₘ))) *
                inv(2.64 - 2.52 * inv(1 + exp(-4.8 - 0.04Vₘ)))

            du[4] =
                0.7462686567164178 * (1.0 + exp(-6.29 - 0.1Vₘ)) *
                (-NaV₊h + inv(1.0 + exp(9.44015444015444 + 0.19305019305019305Vₘ))) *
                inv(1.5 + inv(1.0 + exp(9.694444444444445 + 0.2777777777777778Vₘ)))

            du[5] =
                (-CaT₊m + inv(1.0 + exp(-3.7638888888888893 + -0.1388888888888889Vₘ))) *
                inv(43.4 - 42.6 * inv(1.0 + exp(-3.321951219512195 - 0.04878048780487805Vₘ)))

            du[6] =
                (-CaT₊h + inv(1.0 + exp(5.836363636363637 + 0.18181818181818182Vₘ))) *
                inv(210.0 - 179.6 * inv(1.0 + exp(-3.2544378698224854 - 0.0591715976331361Vₘ)))

            du[7] =
                (-CaS₊m + inv(1.0 + exp(-4.074074074074074 - 0.1234567901234568Vₘ))) *
                inv(2.8 + 14.0 * inv(exp(-5.384615384615385 - 0.07692307692307693Vₘ) + exp(2.7 + 0.1Vₘ)))

            du[8] =
                (-CaS₊h + inv(1.0 + exp(9.67741935483871 + 0.16129032258064516Vₘ))) *
                inv(120.0 + 300.0 * inv(exp(-4.0625 - 0.0625Vₘ) +
                                        exp(6.111111111111111 + 0.1111111111111111Vₘ)))

            du[9] =
                (-KA₊m + inv(1.0 + exp(-3.1264367816091956 - 0.1149425287356322Vₘ))) *
                inv(23.2 - 20.8 * inv(1.0 + exp(-2.164473684210526 - 0.06578947368421052Vₘ)))

            du[10] =
                (-KA₊h + inv(1.0 + exp(11.612244897959181 + 0.2040816326530612Vₘ))) *
                inv(77.2 - 58.4 * inv(1.0 + exp(-1.4679245283018867 - 0.03773584905660377Vₘ)))
            du[11] =
                (-KCa₊m + Caᵢ * inv(3.0 + Caᵢ) *
                inv(1.0 + exp(-2.246031746031746 - 0.07936507936507936Vₘ))) *
                inv(180.6 + -150.2 * inv(1.0 + exp(-2.026431718061674 - 0.04405286343612335Vₘ)))

            du[12] =
                (-Kdr₊m + inv(1.0 + exp(-1.0423728813559323 - 0.0847457627118644Vₘ))) *
                inv(14.4 - 12.8 * inv(1.0 + exp(-1.4739583333333335 - 0.052083333333333336Vₘ)))

            du[13] =
                1 // 2 * 
                (exp(-1.8671328671328669 + 0.06993006993006992Vₘ) +
                 exp(-14.629310344827585 - 0.08620689655172414Vₘ)) *
                (-H₊m + inv(1.0 + exp(13.636363636363637 + 0.18181818181818182Vₘ)))
            nothing
        end
    end
end

byhand_prob = ODEProblem{true}(prinz_neuron!, u0, (0., time), p)
mtk_prob = ODEProblem(simul_sys, [], (0., time), [])
byhand_sol = solve(byhand_prob, Rosenbrock23(), reltol=1e-8, abstol=1e-8);
current_mtk_sol = solve(mtk_prob, Rosenbrock23(), reltol=1e-8, abstol=1e-8);

tsteps = 0.0:0.025:2000.0
byhand_out = Array(byhand_sol(tsteps, idxs=2))
current_mtk_out = current_mtk_sol(tsteps)[Vₘ]

# FIXME: tolerance should need to be this high
@test isapprox(byhand_out, current_mtk_out, rtol=0.1)

end # testset
end # module
