# Example of writing synaptic kinetics
using Conductor, IfElse, OrdinaryDiffEq, Plots, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K

Vₘ = MembranePotential()

nav_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]

kdr_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         0.125 * exp(-(Vₘ + 65.0)/80.0),
         4, name = :n)]

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS/cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS/cm^2)
channels = [NaV, Kdr, leak];
# Equilibrium potentials are a implicit description of a ion concentration gradient
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named Iₑ = IonCurrent(NonIonic)
holding_current = Iₑ ~ ustrip(Float64, µA, 400pA)

@named neuron1 = Compartment(Vₘ, channels, reversals;
                             geometry = Sphere(radius = 20µm),
                             stimuli = [holding_current])

@named neuron2 = Compartment(Vₘ, channels, reversals;
                             geometry = Sphere(radius = 20µm))
                                   
# Synaptic model
syn∞ = 1/(1 + exp((-35 - Vₘ)/5))
τsyn = (1 - syn∞)/(1/40)
syn_kinetics = Gate(SteadyStateTau, syn∞, τsyn, name = :s)
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@named Glut = SynapticChannel(Cation, [syn_kinetics]; max_s = 30nS);

network = NeuronalNetworkSystem([Synapse(neuron1 => neuron2, Glut, EGlut)])

t = 250
sim = Simulation(network, time = t*ms)
solution = solve(sim, Rosenbrock23())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(t*5))

