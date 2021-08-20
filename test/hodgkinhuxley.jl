
module HodgkinHuxleySingle

using Test
using Conductor, ModelingToolkit, OrdinaryDiffEq, Unitful, IfElse
import Conductor: Na, K
using Unitful: mV, mS, cm, µm, µA, ms

@testset "Hodgkin Huxley Single Compartment" begin

# Symbolic model construction

Vₘ = MembranePotential()

nav_kinetics = [
    Gate(AlphaBetaRates,
         αₘ = IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         βₘ = 4.0*exp(-(Vₘ + 65.0)/18.0), p = 3)
    Gate(AlphaBetaRates,
         αₕ = 0.07*exp(-(Vₘ+65.0)/20.0),
         βₕ = 1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)))]

kdr_kinetics = [
    Gate(AlphaBetaRates,
         αₙ = IfElse.ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         βₙ = 0.125 * exp(-(Vₘ + 65.0)/80.0), p = 4)]

@named NaV = IonChannel(Sodium, nav_kinetics, 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, 36mS/cm^2)
@named leak = PassiveChannel(Leak, 0.3mS/cm^2)

gradients = Equilibria([Na   =>  50.0mV, K    => -77.0mV, Leak => -54.4mV])
area = 4*pi*(20µm)^2
t = 300.
pulse(t, current) = IfElse.ifelse(t > 100.0, IfElse.ifelse(t < 200.0, 0.0004, 0.0), 0.0)

@named neuron = Soma([NaV,Kdr,leak], gradients, stimulus = pulse, area = ustrip(Float64, cm^2, area));
sim = Simulation(neuron, time = t*ms) # ODESystem -> structural_simplify -> ODAEProblem output

# TODO: Add test to check that the model contains generated features that we expect

### Hand-written model setup ###
# Steady-state functions
minf(V0) = inv((4.0 + 0.1V0)*inv(1.0 - exp(-4.0 - 0.1V0)) + 4.0*exp(-3.6111111111111107 - 0.05555555555555555V0)) *
           (4.0 + 0.1V0) * inv(1.0 - exp(-4.0 - 0.1V0))
ninf(V0) = (0.55 + 0.01V0) * inv(1.0 - exp(-5.5 - 0.1V0)) * inv(0.125exp(-0.8125 - 0.0125V0) +
                                                                (0.55 + 0.01V0) * inv(1.0 - exp(-5.5 - 0.1V0)))
hinf(V0) = 0.07*exp(-3.25 - 0.05V0) * inv(0.07*exp(-3.25 - 0.05V0) + inv(1.0 + exp(-3.5 - 0.1V0)) )
hand_pulse(t) = ifelse(100. < t < 200. , 0.0004 , 0.0)

# Initial conditions
V0 = -65.0
u0 = zeros(Float64, 4)
p = zeros(Float64, 8)

# Serialized MTK function initial conditions
u0[1] = V0       # Vₘ    
u0[2] = minf(V0) # NaV₊m 
u0[3] = hinf(V0) # NaV₊h 
u0[4] = ninf(V0) # Kdr₊n 

# Parameters 
p[1]  = 0.001      # cₘ      
p[2]  = 5.02655e-5 # aₘ      
p[3]  = 50.0       # ENa     
p[4]  = -77.0      # EK      
p[5]  = -54.4      # El      
p[6]  = 120.0      # NaV₊gbar
p[7]  = 36.0       # Kdr₊gbar
p[8]  = 0.3        # leak₊g  

# Cleaned up version of what default MTK output looks like
function hodgkin_huxley!(du, u, p, t)
    @inbounds let Vₘ         = u[1],
                  NaV₊m      = u[2],
                  NaV₊h      = u[3],
                  Kdr₊n      = u[4],
          
                  cₘ          = p[1],
                  aₘ          = p[2],
                  ENa         = p[3],
                  EK          = p[4],
                  El          = p[5],
                  NaV₊gbar    = p[6],
                  Kdr₊gbar    = p[7],
                  leak₊g      = p[8],
                  Iapp        = hand_pulse(t)
    du[1] = (Iapp - aₘ*leak₊g*(Vₘ - El) - Kdr₊gbar*aₘ*Kdr₊n^4*(Vₘ - EK) -
             NaV₊gbar*aₘ*NaV₊h*NaV₊m^3*(Vₘ - ENa))*inv(aₘ)*inv(cₘ)
    du[2] = ifelse(Vₘ == -40.0 , 1.0 , (4.0 + 0.1Vₘ) * inv(1.0 - exp(-4.0 - 0.1Vₘ))) * 
            (1 - NaV₊m) - 4.0*NaV₊m*exp(-3.6111111111111107 + -0.05555555555555555Vₘ)

    du[3] = 0.07*exp(-3.25 - 0.05Vₘ)*(1 - NaV₊h) - NaV₊h*inv(1.0 + exp(-3.5 - 0.1Vₘ))
    du[4] = ifelse(Vₘ == -55.0 , 0.1 , (0.55 + 0.01Vₘ) * inv(1.0 - exp(-5.5 - 0.1Vₘ)))*(1 - Kdr₊n) - 
            0.125 * Kdr₊n * exp(-0.8125 - 0.0125Vₘ)
    end
    return nothing
end

byhand_prob = ODEProblem{true}(hodgkin_huxley!, u0, (0.,300.), p)
byhand_sol = solve(byhand_prob, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);
current_mtk_sol = solve(sim, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);

tsteps = 0.0:0.025:300.0
byhand_out = Array(byhand_sol(tsteps, idxs=1))
current_mtk_out = current_mtk_sol(tsteps)[Vₘ]

@test isapprox(byhand_out, current_mtk_out, rtol=0.0001)

end #testset

end # module
