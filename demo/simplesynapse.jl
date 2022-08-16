# Example of writing synaptic kinetics
using Conductor, IfElse, OrdinaryDiffEq, Plots, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K

include("hh_rates.jl")
Vₘ = MembranePotential(-65mV)

nav_kinetics = [Gate(AlphaBeta, αₘ(Vₘ), βₘ(Vₘ), p = 3, name = :m),
                Gate(AlphaBeta, αₕ(Vₘ), βₕ(Vₘ), name = :h)]
kdr_kinetics = [Gate(AlphaBeta, αₙ(Vₘ), βₙ(Vₘ), p = 4, name = :n)]

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named Iₑ = IonCurrent(NonIonic)
@named I_hold = IonCurrent(NonIonic, 5000pA, dynamic = false)
holding_current = Iₑ ~ I_hold

dynamics_1 = HodgkinHuxley(Vₘ, channels, reversals;
                           geometry = Cylinder(radius = 25µm, height = 400µm),
                           stimuli = [holding_current]);
dynamics_2 = HodgkinHuxley(Vₘ, channels, reversals;
                           geometry = Cylinder(radius = 25µm, height = 400µm));

@named neuron1 = Compartment(dynamics_1)
@named neuron2 = Compartment(dynamics_2)
                                   
# Synaptic model
Vₓ = ExtrinsicPotential()

syn_kinetics = Gate(SteadyStateTau, syn∞(Vₓ), τsyn(Vₓ), name = :z)
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@named Glut = SynapticChannel(Cation, [syn_kinetics]; max_s = 30nS);

topology = NetworkTopology([neuron1, neuron2], [Glut]);
topology[neuron1, neuron2] = Glut
reversal_map = Dict([Glut => EGlut])

@named net = NeuronalNetworkSystem(topology, reversal_map)
total_time = 250.0
sim = Simulation(net, time = total_time*ms)

solution = solve(sim, Rosenbrock23(), abstol=1e-3, reltol=1e-3, saveat=0.2)

# Plot at 5kHz sampling
plot(solution(0.0:0.2:total_time, idxs = [neuron1.Vₘ, neuron2.Vₘ]))

