
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__); Pkg.instantiate()

include("prinz_neuron_setup.jl")
include("helpers.jl")

system = neuron.sys
@named simulation = ODESystem([D(system.Isyn) ~ 0]; systems = [system])
simplified = structural_simplify(simulation)

sim_exp = ModelingToolkit.build_torn_function(simplified; expression=true)
cleaned = clean_expr!(sim_exp, simplified)

clipboard(format_text(string(prettify(cleaned)))) # can paste readable output below

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
                  EK          = p[8],
                  EK          = p[9],
                  EH          = p[10],
                  El          = p[11],
                  NaV₊gbar    = p[12],
                  CaT₊gbar    = p[13],
                  CaS₊gbar    = p[14],
                  KA₊gbar     = p[15],
                  KCa₊gbar    = p[16],
                  Kdr₊gbar    = p[17],
                  H₊gbar      = p[18],
                  leak₊g      = p[19]

        let
            du[1] =
                (
                    Ca∞ +
                    -1Caᵢ +
                    -1000 *
                    f *
                    (
                        CaS₊gbar *
                        aₘ *
                        CaS₊h *
                        CaS₊m^3.0 *
                        (
                            Vₘ +
                            12.199306599503455 *
                            log(0.0003333333333333333 * max(Caᵢ, 0.001))
                        ) +
                        CaT₊gbar *
                        aₘ *
                        CaT₊h *
                        CaT₊m^3.0 *
                        (
                            Vₘ +
                            12.199306599503455 *
                            log(0.0003333333333333333 * max(Caᵢ, 0.001))
                        )
                    )
                ) * inv(τCa)
            du[2] =
                (
                    -1 * aₘ * leak₊g * (Vₘ + -1El) +
                    -1 *
                    H₊gbar *
                    aₘ *
                    H₊m *
                    (Vₘ + -1EH) +
                    -1 *
                    KCa₊gbar *
                    aₘ *
                    KCa₊m^4.0 *
                    (Vₘ + -1EK) +
                    -1 *
                    Kdr₊gbar *
                    aₘ *
                    Kdr₊m^4.0 *
                    (Vₘ + -1EK) +
                    -1 *
                    CaS₊gbar *
                    aₘ *
                    CaS₊h *
                    CaS₊m^3.0 *
                    (
                        Vₘ +
                        12.199306599503455 *
                        log(0.0003333333333333333 * max(Caᵢ, 0.001))
                    ) +
                    -1 *
                    CaT₊gbar *
                    aₘ *
                    CaT₊h *
                    CaT₊m^3.0 *
                    (
                        Vₘ +
                        12.199306599503455 *
                        log(0.0003333333333333333 * max(Caᵢ, 0.001))
                    ) +
                    -1 *
                    KA₊gbar *
                    aₘ *
                    KA₊h *
                    KA₊m^3.0 *
                    (Vₘ + -1EK) +
                    -1 *
                    NaV₊gbar *
                    aₘ *
                    NaV₊h *
                    NaV₊m^3.0 *
                    (Vₘ + -1ENa)
                ) *
                inv(aₘ) *
                inv(cₘ)
            du[3] =
                (
                    -1NaV₊m +
                    inv(1.0 + exp(-4.8204158790170135 + -0.1890359168241966Vₘ))
                ) * inv(2.64 + -2.52 * inv(1 + exp(-4.8 + -0.04Vₘ)))
            du[4] =
                0.7462686567164178 *
                (1.0 + exp(-6.29 + -0.1Vₘ)) *
                (
                    -1NaV₊h +
                    inv(1.0 + exp(9.44015444015444 + 0.19305019305019305Vₘ))
                ) *
                inv(1.5 + inv(1.0 + exp(9.694444444444445 + 0.2777777777777778Vₘ)))
            du[5] =
                (
                    -1CaT₊m +
                    inv(1.0 + exp(-3.7638888888888893 + -0.1388888888888889Vₘ))
                ) * inv(
                    43.4 +
                    -42.6 *
                    inv(1.0 + exp(-3.321951219512195 + -0.04878048780487805Vₘ)),
                )
            du[6] =
                (
                    -1CaT₊h +
                    inv(1.0 + exp(5.836363636363637 + 0.18181818181818182Vₘ))
                ) * inv(
                    210.0 +
                    -179.6 *
                    inv(1.0 + exp(-3.2544378698224854 + -0.0591715976331361Vₘ)),
                )
            du[7] =
                (
                    -1CaS₊m +
                    inv(1.0 + exp(-4.074074074074074 + -0.1234567901234568Vₘ))
                ) * inv(
                    2.8 +
                    14.0 * inv(
                        exp(-5.384615384615385 + -0.07692307692307693Vₘ) +
                        exp(2.7 + 0.1Vₘ),
                    ),
                )
            du[8] =
                (
                    -1CaS₊h +
                    inv(1.0 + exp(9.67741935483871 + 0.16129032258064516Vₘ))
                ) * inv(
                    120.0 +
                    300.0 * inv(
                        exp(-4.0625 + -0.0625Vₘ) +
                        exp(6.111111111111111 + 0.1111111111111111Vₘ),
                    ),
                )
            du[9] =
                (
                    -1KA₊m +
                    inv(1.0 + exp(-3.1264367816091956 + -0.1149425287356322Vₘ))
                ) * inv(
                    23.2 +
                    -20.8 *
                    inv(1.0 + exp(-2.164473684210526 + -0.06578947368421052Vₘ)),
                )
            du[10] =
                (
                    -1KA₊h +
                    inv(1.0 + exp(11.612244897959181 + 0.2040816326530612Vₘ))
                ) * inv(
                    77.2 +
                    -58.4 *
                    inv(1.0 + exp(-1.4679245283018867 + -0.03773584905660377Vₘ)),
                )
            du[11] =
                (
                    -1KCa₊m +
                    Caᵢ *
                    inv(3.0 + Caᵢ) *
                    inv(1.0 + exp(-2.246031746031746 + -0.07936507936507936Vₘ))
                ) * inv(
                    180.6 +
                    -150.2 *
                    inv(1.0 + exp(-2.026431718061674 + -0.04405286343612335Vₘ)),
                )
            du[12] =
                (
                    -1Kdr₊m +
                    inv(1.0 + exp(-1.0423728813559323 + -0.0847457627118644Vₘ))
                ) * inv(
                    14.4 +
                    -12.8 *
                    inv(1.0 + exp(-1.4739583333333335 + -0.052083333333333336Vₘ)),
                )
            du[13] =
                1 // 2 *
                (
                    exp(-1.8671328671328669 + 0.06993006993006992Vₘ) +
                    exp(-14.629310344827585 + -0.08620689655172414Vₘ)
                ) *
                (
                    -1H₊m +
                    inv(1.0 + exp(13.636363636363637 + 0.18181818181818182Vₘ))
                )
            nothing
        end
    end
end
