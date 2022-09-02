
module HodgkinHuxleySingle

using Test
using Conductor, ModelingToolkit, OrdinaryDiffEq, Unitful
import ModelingToolkit: isparameter
import Conductor: Na, K
using Unitful: mV, mS, cm, µm, µA, ms, pA

@testset "Hodgkin Huxley Single Compartment" begin

# Symbolic model construction

Vₘ = MembranePotential(-65mV)

nav_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

kdr_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4, name = :n)]

alphamss = ((4.0exp(-3.6111111111111107 - (0.05555555555555555Vₘ)) +
            ifelse(Vₘ == -40.0, 1.0, (4.0 + 0.1Vₘ)*((1.0 - exp(-4.0 - (0.1Vₘ)))^-1)))^-1) * 
            ifelse(Vₘ == -40.0, 1.0, (4.0 + 0.1Vₘ)*((1.0 - exp(-4.0 - (0.1Vₘ)))^-1))

@test isequal(steadystate(nav_kinetics[1]), alphamss)

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

@test isparameter(leak.gbar) && !isparameter(get_output(NaV))
channels = [NaV, Kdr, leak]
@test [length(equations(x)) for x in channels] == [3,2,1]
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])
@named Iₑ = IonCurrent(NonIonic)
@named I_rest = IonCurrent(NonIonic, 0.0µA, dynamic = false)
@named I_step = IonCurrent(NonIonic, 400.0pA, dynamic = false)
@parameters tstart = 100.0 [unit=ms] tstop = 200.0 [unit=ms]
electrode_pulse = Iₑ ~ ifelse((t > tstart) & (t < tstop), I_step, I_rest)

dynamics = HodgkinHuxley(Vₘ, channels, reversals;
                         geometry = Sphere(radius = 20µm),
                         stimuli = [electrode_pulse])

@named neuron = Compartment(dynamics)

@test length.([equations(neuron),
               states(neuron),
               parameters(neuron)]) == [13,13,12]

time = 300.
sim_sys = Simulation(neuron, time = time*ms, return_system = true)

@test length.([equations(sim_sys),
               states(sim_sys),
               parameters(sim_sys)]) == [4,4,12]

expect_syms = [:Vₘ, :NaV₊m, :NaV₊h, :Kdr₊n]
@test all(x -> hasproperty(sim_sys, x), expect_syms)

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
p[1]  = 1.0        # cₘ      
#p[2]  = 5.02655e-5 # aₘ      
p[2]  = area(Sphere(radius = 20µm)) # aₘ      
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
mtk_prob = ODEProblem(sim_sys, [], (0., time), [])

byhand_sol = solve(byhand_prob, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);
current_mtk_sol = solve(mtk_prob, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);

tsteps = 0.0:0.025:300.0
byhand_out = Array(byhand_sol(tsteps, idxs=1))
current_mtk_out = current_mtk_sol(tsteps)[Vₘ]

@test isapprox(byhand_out, current_mtk_out, rtol=0.0001)

end #testset

end # module
