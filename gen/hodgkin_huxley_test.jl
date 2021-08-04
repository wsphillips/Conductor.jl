
V0 = -65.0
u0 = zeros(Float64, 6)
u  = zeros(Float64, 6)
p  = zeros(Float64, 8)
du = similar(u)
# Steady-state functions
minf(V0) = inv((4.0 + 0.1V0)*inv(1.0 - exp(-4.0 - 0.1V0)) + 4.0*exp(-3.6111111111111107 - 0.05555555555555555V0)) * (4.0 + 0.1V0)*inv(1.0 - exp(-4.0 - 0.1V0))
ninf(V0) = (0.55 + 0.01V0) * inv(1.0 - exp(-5.5 - 0.1V0)) * inv(0.125exp(-0.8125 - 0.0125V0) + (0.55 + 0.01V0) * inv(1.0 - exp(-5.5 - 0.1V0)))
hinf(V0) = 0.07*exp(-3.25 - 0.05V0) * inv(0.07*exp(-3.25 - 0.05V0) + inv(1.0 + exp(-3.5 - 0.1V0)) )
pulse(t, current) = 100. < t < 200. ? 0.0004 : 0.0

# Initial conditions
u0[1] = 0        # Isyn 
u0[2] = V0       # Vₘ   
u0[3] = minf(V0) # NaV₊m
u0[4] = hinf(V0) # NaV₊h
u0[5] = ninf(V0) # Kdr₊n
u0[6] = 0.0      # Iapp 

# Parameters 
p[1]  = 0.001      # cₘ      
p[2]  = 5.02655e-5 # aₘ      
p[3]  = 50.0       # ENa     
p[4]  = -77.0      # EK      
p[5]  = -54.4      # El      
p[6]  = 120.0      # NaV₊gbar
p[7]  = 36.0       # Kdr₊gbar
p[8]  = 0.3        # leak₊g  

# Layout
diff_vars = [true, true, true, true, true, false]

function hodgkin_huxley!(du, u, p, t)

    @inbounds let Isyn        = u[1],
                  Vₘ          = u[2],
                  NaV₊m       = u[3],
                  NaV₊h       = u[4],
                  Kdr₊n       = u[5],
                  Iapp        = u[6],
          
                  cₘ          = p[1],
                  aₘ          = p[2],
                  ENa         = p[3],
                  EK          = p[4],
                  El          = p[5],
                  NaV₊gbar    = p[6],
                  Kdr₊gbar    = p[7],
                  leak₊g      = p[8]

    du[1] = 0
    du[2] = (Iapp - aₘ*leak₊g*(Vₘ - El) - Kdr₊gbar*aₘ*Kdr₊n^4.0*(Vₘ - EK) - Isyn -
             NaV₊gbar*aₘ*NaV₊h*NaV₊m^3.0*(Vₘ - ENa))*inv(aₘ)*inv(cₘ)

    du[3] = (Vₘ == -40.0 ? 1.0 : (4.0 + 0.1Vₘ) * inv(1.0 - exp(-4.0 + -0.1Vₘ))) * 
            (1 - NaV₊m) - 4.0*NaV₊m*exp(-3.6111111111111107 + -0.05555555555555555Vₘ)

    du[4] = 0.07*exp(-3.25 + -0.05Vₘ)*(1 - NaV₊h) - NaV₊h*inv(1.0 + exp(-3.5 + -0.1Vₘ))

    du[5] = (Vₘ == -55.0 ? 0.1 : (0.55 + 0.01Vₘ) * inv(1.0 - exp(-5.5 + -0.1Vₘ)))*(1 - Kdr₊n) - 
            0.125 * Kdr₊n * exp(-0.8125 - 0.0125Vₘ)

    du[6] = pulse(t, Iapp) - Iapp
    end
    return nothing
end

function raw_hodgkin_huxley!(du, u, p, t)
    let neuron₊Isyn = @inbounds(u[1]),
        neuron₊Vₘ = @inbounds(u[2]),
        neuron₊NaV₊m = @inbounds(u[3]),
        neuron₊NaV₊h = @inbounds(u[4]),
        neuron₊Kdr₊n = @inbounds(u[5]),
        neuron₊Iapp = @inbounds(u[6]),
        neuron₊cₘ = @inbounds(p[1]),
        neuron₊aₘ = @inbounds(p[2]),
        neuron₊ENa = @inbounds(p[3]),
        neuron₊EK = @inbounds(p[4]),
        neuron₊El = @inbounds(p[5]),
        neuron₊NaV₊gbar = @inbounds(p[6]),
        neuron₊Kdr₊gbar = @inbounds(p[7]),
        neuron₊leak₊g = @inbounds(p[8])

        @inbounds begin
            du[1] = 0
            du[2] = pulse(t, neuron₊Iapp) + -1neuron₊Iapp
            du[3] =
                (
                    (
                        -1 * neuron₊aₘ * neuron₊leak₊g * (neuron₊Vₘ + -1neuron₊El) +
                        (-1 * neuron₊Kdr₊gbar * neuron₊aₘ * neuron₊Kdr₊n^4.0) *
                        (*)(neuron₊Vₘ + -1neuron₊EK) +
                        neuron₊Iapp +
                        -1neuron₊Isyn
                    ) + +(
                        (-1 * neuron₊NaV₊gbar * neuron₊aₘ * neuron₊NaV₊h) *
                        (neuron₊NaV₊m^3.0 * (neuron₊Vₘ + -1neuron₊ENa)),
                    )
                ) *
                inv(neuron₊aₘ) *
                inv(neuron₊cₘ)
            du[4] =
                if neuron₊Vₘ == -40.0
                    1.0
                else
                    (4.0 + 0.1neuron₊Vₘ) * inv(1.0 + -1 * exp(-4.0 + -0.1neuron₊Vₘ))
                end * (1 + -1neuron₊NaV₊m) +
                -4.0 *
                neuron₊NaV₊m *
                exp(-3.6111111111111107 + -0.05555555555555555neuron₊Vₘ)
            du[5] =
                -1 * neuron₊NaV₊h * inv(1.0 + exp(-3.5 + -0.1neuron₊Vₘ)) +
                0.07 * exp(-3.25 + -0.05neuron₊Vₘ) * (1 + -1neuron₊NaV₊h)
            du[6] =
                (1 + -1neuron₊Kdr₊n) * if neuron₊Vₘ == -55.0
                    0.1
                else
                    (0.55 + 0.01neuron₊Vₘ) * inv(1.0 + -1 * exp(-5.5 + -0.1neuron₊Vₘ))
                end + -0.125 * neuron₊Kdr₊n * exp(-0.8125 + -0.0125neuron₊Vₘ)
            nothing
        end
    end
end
#=
function hodgkin_huxley!(du, u, p, t)

    let Isyn        = @inbounds(u[1]),
        Vₘ          = @inbounds(u[2]),
        NaV₊m       = @inbounds(u[3]),
        NaV₊h       = @inbounds(u[4]),
        Kdr₊n       = @inbounds(u[5]),
        Iapp        = @inbounds(u[6]),

        cₘ          = @inbounds(p[1]),
        aₘ          = @inbounds(p[2]),
        ENa         = @inbounds(p[3]),
        EK          = @inbounds(p[4]),
        El          = @inbounds(p[5]),
        NaV₊gbar    = @inbounds(p[6]),
        Kdr₊gbar    = @inbounds(p[7]),
        leak₊g      = @inbounds(p[8])

        @inbounds begin
            du[1] = 0
            du[2] = pulse(t, Iapp) - Iapp

            du[3] = (Iapp - aₘ*leak₊g*(Vₘ - El) - Kdr₊gbar*aₘ*Kdr₊n^4.0*(Vₘ - EK) - Isyn -
                     NaV₊gbar*aₘ*NaV₊h*NaV₊m^3.0*(Vₘ - ENa))*inv(aₘ)*inv(cₘ)

            du[4] = (Vₘ == -40.0 ? 1.0 : (4.0 + 0.1Vₘ) * inv(1.0 - exp(-4.0 + -0.1Vₘ))) * 
                    (1 - NaV₊m) - 4.0*NaV₊m*exp(-3.6111111111111107 + -0.05555555555555555Vₘ)

            du[5] = 0.07*exp(-3.25 + -0.05Vₘ)*(1 - NaV₊h) - NaV₊h*inv(1.0 + exp(-3.5 + -0.1Vₘ))

            du[6] = (Vₘ == -55.0 ? 0.1 : (0.55 + 0.01Vₘ) * inv(1.0 - exp(-5.5 + -0.1Vₘ)))*(1 - Kdr₊n) - 
                    0.125 * Kdr₊n * exp(-0.8125 - 0.0125Vₘ)
            nothing
        end
    end
end

=#
