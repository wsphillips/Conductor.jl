
module SimpleSynapse

using Test, Conductor, OrdinaryDiffEq, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K

@testset "Simple excitatory synapse" begin

Vₘ = MembranePotential(-65mV)

nav_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0),
         p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

kdr_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4, name = :n)]

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named Iₑ = IonCurrent(NonIonic)
@named I_hold = IonCurrent(NonIonic, 5000pA, dynamic = false)
holding_current = Iₑ ~ I_hold

geo = Cylinder(radius = 25µm, height = 400µm)

dynamics_1 = HodgkinHuxley(Vₘ, channels, reversals; geometry = geo, stimuli = [holding_current]);
dynamics_2 = HodgkinHuxley(Vₘ, channels, reversals; geometry = geo);
@named neuron1 = Compartment(dynamics_1)
@named neuron2 = Compartment(dynamics_2)
 
# Synaptic model
Vₓ = ExtrinsicPotential()
syn∞ = 1/(1 + exp((-35 - Vₓ)/5))
τsyn = (1 - syn∞)/(1/40)
syn_kinetics = Gate(SteadyStateTau, syn∞, τsyn, name = :z)
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@named Glut = SynapticChannel(Cation, [syn_kinetics]; max_s = 30nS);

@test length.([equations(Glut),
               states(Glut),
               parameters(Glut)]) == [2,3,1]

topology = NetworkTopology([neuron1, neuron2], [Glut]);
topology[neuron1, neuron2] = Glut
reversal_map = Dict([Glut => EGlut])

@named network = NeuronalNetworkSystem(topology, reversal_map)

@test length.([equations(network),
               states(network),
               parameters(network)]) == [29,29,19]

ttot = 250.
simul_sys = Simulation(network, time = ttot*ms, return_system = true)

@test length.([equations(simul_sys),
               states(simul_sys),
               parameters(simul_sys)]) == [9,9,19]

# Skipped because synaptic systems generated via Base.gensym
#expect_syms = [:neuron1₊Isyn] # just take one for now
#@test all(x -> hasproperty(simul_sys, x), expect_syms)

# Hand written problem
# Steady-state functions
minf(V0) = inv((4.0 + 0.1V0)*inv(1.0 - exp(-4.0 - 0.1V0)) + 4.0*exp(-3.6111111111111107 - 0.05555555555555555V0)) *
           (4.0 + 0.1V0) * inv(1.0 - exp(-4.0 - 0.1V0))
ninf(V0) = (0.55 + 0.01V0) * inv(1.0 - exp(-5.5 - 0.1V0)) * inv(0.125exp(-0.8125 - 0.0125V0) +
                                                                (0.55 + 0.01V0) * inv(1.0 - exp(-5.5 - 0.1V0)))
hinf(V0) = 0.07*exp(-3.25 - 0.05V0) * inv(0.07*exp(-3.25 - 0.05V0) + inv(1.0 + exp(-3.5 - 0.1V0)) )

# Initial conditions
V0 = -65.0
u0 = zeros(Float64, 12)
p = zeros(Float64, 18)

# Serialized MTK function initial conditions
# Neuron 1
u0[1]  = 0        # neuron1₊Isyn 
u0[2]  = 0.005    # neuron1₊Iapp 
u0[3]  = V0       # neuron1₊Vₘ   
u0[4]  = minf(V0) # neuron1₊NaV₊m
u0[5]  = hinf(V0) # neuron1₊NaV₊h
u0[6]  = ninf(V0) # neuron1₊Kdr₊n

# Neuron 2
u0[7]  = 0.0      # neuron2₊Iapp 
u0[8]  = V0       # neuron2₊Vₘ   
u0[9]  = minf(V0) # neuron2₊NaV₊m
u0[10] = hinf(V0) # neuron2₊NaV₊h
u0[11] = ninf(V0) # neuron2₊Kdr₊n

# Synapse
u0[12] = 0.0      # Glut1₊s      

# Parameters 
# Neuron 1
p[2]   =  1.0       # neuron1₊cₘ      
p[3]   =  area(geo) # neuron1₊aₘ      
p[4]   =  50.0      # neuron1₊ENa     
p[5]   =  -77.0     # neuron1₊EK      
p[6]   =  -54.4     # neuron1₊El      
p[7]   =  120.0     # neuron1₊NaV₊gbar
p[8]   =  36.0      # neuron1₊Kdr₊gbar
p[9]   =  0.3       # neuron1₊leak₊g  

# Neuron 2
p[10]  =  1.0       # neuron2₊cₘ      
p[11]  =  0.000629  # neuron2₊aₘ      
p[12]  =  50.0      # neuron2₊ENa     
p[13]  =  -77.0     # neuron2₊EK      
p[14]  =  -54.4     # neuron2₊El      
p[15]  =  120.0     # neuron2₊NaV₊gbar
p[16]  =  36.0      # neuron2₊Kdr₊gbar
p[17]  =  0.3       # neuron2₊leak₊g  

# Synapse
p[1]   =  0.0       # EGlut           
p[18]  =  3.0e-5    # Glut1₊gbar      

function simple_synapse!(du, u, p, t)
    let neuron1₊Isyn        = @inbounds(u[1]),
        neuron1₊Iapp        = @inbounds(u[2]),
        neuron1₊Vₘ          = @inbounds(u[3]),
        neuron1₊NaV₊m       = @inbounds(u[4]),
        neuron1₊NaV₊h       = @inbounds(u[5]),
        neuron1₊Kdr₊n       = @inbounds(u[6]),
        neuron2₊Iapp        = @inbounds(u[7]),
        neuron2₊Vₘ          = @inbounds(u[8]),
        neuron2₊NaV₊m       = @inbounds(u[9]),
        neuron2₊NaV₊h       = @inbounds(u[10]),
        neuron2₊Kdr₊n       = @inbounds(u[11]),
        Glut1₊s             = @inbounds(u[12]),

        EGlut               = @inbounds(p[1]),
        neuron1₊cₘ          = @inbounds(p[2]),
        neuron1₊aₘ          = @inbounds(p[3]),
        neuron1₊ENa         = @inbounds(p[4]),
        neuron1₊EK          = @inbounds(p[5]),
        neuron1₊El          = @inbounds(p[6]),
        neuron1₊NaV₊gbar    = @inbounds(p[7]),
        neuron1₊Kdr₊gbar    = @inbounds(p[8]),
        neuron1₊leak₊g      = @inbounds(p[9]),
        neuron2₊cₘ          = @inbounds(p[10]),
        neuron2₊aₘ          = @inbounds(p[11]),
        neuron2₊ENa         = @inbounds(p[12]),
        neuron2₊EK          = @inbounds(p[13]),
        neuron2₊El          = @inbounds(p[14]),
        neuron2₊NaV₊gbar    = @inbounds(p[15]),
        neuron2₊Kdr₊gbar    = @inbounds(p[16]),
        neuron2₊leak₊g      = @inbounds(p[17]),
        Glut1₊gbar          = @inbounds(p[18])

        let
            du[1] = 0
            du[2] = 0
            du[3] = (-1 * neuron1₊aₘ * neuron1₊leak₊g * (neuron1₊Vₘ + -1neuron1₊El) +
                     -1 * neuron1₊Kdr₊gbar * neuron1₊aₘ * neuron1₊Kdr₊n^4.0 * (neuron1₊Vₘ + -1neuron1₊EK) +
                     -1 * neuron1₊NaV₊gbar * neuron1₊aₘ * neuron1₊NaV₊h * neuron1₊NaV₊m^3.0 * (neuron1₊Vₘ + -1neuron1₊ENa) +
                     neuron1₊Iapp + -1neuron1₊Isyn) * inv(neuron1₊aₘ) * inv(neuron1₊cₘ)

            du[4] = if neuron1₊Vₘ == -40.0
                1.0
            else
                (4.0 + 0.1neuron1₊Vₘ) * inv(1.0 + -1 * exp(-4.0 + -0.1neuron1₊Vₘ))
            end * (1 + -1neuron1₊NaV₊m) + -4.0 * neuron1₊NaV₊m *
                    exp(-3.6111111111111107 + -0.05555555555555555neuron1₊Vₘ)

            du[5] = -1 * neuron1₊NaV₊h * inv(1.0 + exp(-3.5 + -0.1neuron1₊Vₘ)) +
                    0.07 * exp(-3.25 + -0.05neuron1₊Vₘ) * (1 + -1neuron1₊NaV₊h)

            du[6] = (1 + -1neuron1₊Kdr₊n) * if neuron1₊Vₘ == -55.0
                0.1
            else
                (0.55 + 0.01neuron1₊Vₘ) * inv(1.0 + -1 * exp(-5.5 + -0.1neuron1₊Vₘ))
            end + -0.125 * neuron1₊Kdr₊n * exp(-0.8125 + -0.0125neuron1₊Vₘ)

            du[7] = 0
            du[8] = (-1 * Glut1₊gbar * Glut1₊s * (neuron2₊Vₘ + -1 * EGlut) +
                     -1 * neuron2₊aₘ * neuron2₊leak₊g * (neuron2₊Vₘ + -1neuron2₊El) +
                     -1 * neuron2₊Kdr₊gbar * neuron2₊aₘ * neuron2₊Kdr₊n^4.0 * (neuron2₊Vₘ + -1neuron2₊EK) +
                     -1 * neuron2₊NaV₊gbar * neuron2₊aₘ * neuron2₊NaV₊h * neuron2₊NaV₊m^3.0 * (neuron2₊Vₘ + -1neuron2₊ENa) +
                     neuron2₊Iapp) * inv(neuron2₊aₘ) * inv(neuron2₊cₘ)

            du[9] = if neuron2₊Vₘ == -40.0
                1.0
            else
                (4.0 + 0.1neuron2₊Vₘ) * inv(1.0 + -1 * exp(-4.0 + -0.1neuron2₊Vₘ))
            end * (1 + -1neuron2₊NaV₊m) + -4.0 * neuron2₊NaV₊m *
                    exp(-3.6111111111111107 + -0.05555555555555555neuron2₊Vₘ)

            du[10] = -1 * neuron2₊NaV₊h * inv(1.0 + exp(-3.5 + -0.1neuron2₊Vₘ)) +
                     0.07 * exp(-3.25 + -0.05neuron2₊Vₘ) * (1 + -1neuron2₊NaV₊h)

            du[11] = (1 + -1neuron2₊Kdr₊n) * if neuron2₊Vₘ == -55.0
                0.1
            else
                (0.55 + 0.01neuron2₊Vₘ) * inv(1.0 + -1 * exp(-5.5 + -0.1neuron2₊Vₘ))
            end + -0.125 * neuron2₊Kdr₊n * exp(-0.8125 + -0.0125neuron2₊Vₘ)

            du[12] = (-1Glut1₊s + inv(1 + exp(-7 // 1 + -1 // 5 * neuron1₊Vₘ))) *
                     inv(40.0 + -40.0 * inv(1 + exp(-7 // 1 + -1 // 5 * neuron1₊Vₘ)))
            nothing
        end
    end
end

# Solve and check for invariance
byhand_prob = ODEProblem{true}(simple_synapse!, u0, (0.,ttot), p)
mtk_prob = ODEProblem(simul_sys, [], (0., ttot), [])
byhand_sol = solve(byhand_prob, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);
current_mtk_sol = solve(mtk_prob, Rosenbrock23(), reltol=1e-9, abstol=1e-9, saveat=0.025);

tsteps = 0.0:0.025:ttot
byhand_out = Array(byhand_sol(tsteps, idxs=8))
current_mtk_out = current_mtk_sol(tsteps)[neuron2.Vₘ]

@test isapprox(byhand_out, current_mtk_out, rtol=0.001)

end # testset

end # module
