
cd(@__DIR__)
using Pkg
Pkg.activate(@__DIR__); Pkg.instantiate()

include("hodgkin_huxley_setup.jl")
include("helpers.jl")

system = neuron.sys
@named simulation = ODESystem([D(system.Isyn) ~ 0]; systems = [system])
simplified = structural_simplify(simulation)

sim_exp = ModelingToolkit.build_torn_function(simplified; expression=true)
cleaned = clean_expr(sim_exp, simplified)
#display_expr(cleaned)
#clipboard(format_text(string(prettify(cleaned)))) # can paste readable output below

############################################################################################
# Generated function expression from MTK lightly cleaned up for better readability

function hh_generated!(du, u, p, t)
    let neuron₊Isyn = @inbounds(u[1]),
        neuron₊Vₘ = @inbounds(u[2]),
        neuron₊NaV₊m = @inbounds(u[3]),
        neuron₊NaV₊h = @inbounds(u[4]),
        neuron₊Kdr₊n = @inbounds(u[5]),
        neuron₊cₘ = @inbounds(p[1]),
        neuron₊aₘ = @inbounds(p[2]),
        neuron₊ENa = @inbounds(p[3]),
        neuron₊EK = @inbounds(p[4]),
        neuron₊El = @inbounds(p[5]),
        neuron₊NaV₊gbar = @inbounds(p[6]),
        neuron₊Kdr₊gbar = @inbounds(p[7]),
        neuron₊leak₊g = @inbounds(p[8])

        let
            du[1] = 0
            du[2] =
                (
                    if t > 100.0
                        if t < 200.0
                            0.0004
                        else
                            0.0
                        end
                    else
                        0.0
                    end +
                    -1neuron₊Isyn +
                    -1 * neuron₊aₘ * neuron₊leak₊g * (neuron₊Vₘ + -1neuron₊El) +
                    -1 *
                    neuron₊Kdr₊gbar *
                    neuron₊aₘ *
                    neuron₊Kdr₊n^4.0 *
                    (neuron₊Vₘ + -1neuron₊EK) +
                    -1 *
                    neuron₊NaV₊gbar *
                    neuron₊aₘ *
                    neuron₊NaV₊h *
                    neuron₊NaV₊m^3.0 *
                    (neuron₊Vₘ + -1neuron₊ENa)
                ) *
                inv(neuron₊aₘ) *
                inv(neuron₊cₘ)
            du[3] =
                (1 + -1neuron₊NaV₊m) * if neuron₊Vₘ == -40.0
                    1.0
                else
                    (4.0 + 0.1neuron₊Vₘ) * inv(1.0 + -1 * exp(-4.0 + -0.1neuron₊Vₘ))
                end +
                -4.0 *
                neuron₊NaV₊m *
                exp(-3.6111111111111107 + -0.05555555555555555neuron₊Vₘ)
            du[4] =
                -1 * neuron₊NaV₊h * inv(1.0 + exp(-3.5 + -0.1neuron₊Vₘ)) +
                0.07 * exp(-3.25 + -0.05neuron₊Vₘ) * (1 + -1neuron₊NaV₊h)
            du[5] =
                if neuron₊Vₘ == -55.0
                    0.1
                else
                    (0.55 + 0.01neuron₊Vₘ) * inv(1.0 + -1 * exp(-5.5 + -0.1neuron₊Vₘ))
                end * (1 + -1neuron₊Kdr₊n) +
                -0.125 * neuron₊Kdr₊n * exp(-0.8125 + -0.0125neuron₊Vₘ)
            nothing
        end
    end
end

