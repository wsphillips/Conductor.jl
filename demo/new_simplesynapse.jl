# Example of writing synaptic kinetics
using Conductor, OrdinaryDiffEq, Plots, Unitful, ModelingToolkit
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K

Vₘ = ParentScope(MembranePotential(-65mV))

nav_kinetics = [Gate(AlphaBeta,
                     ifelse(Vₘ == -40.0, 1.0,
                            (0.1 * (Vₘ + 40.0)) / (1.0 - exp(-(Vₘ + 40.0) / 10.0))),
                     4.0 * exp(-(Vₘ + 65.0) / 18.0),
                     p = 3, name = :m)
                Gate(AlphaBeta,
                     0.07 * exp(-(Vₘ + 65.0) / 20.0),
                     1.0 / (1.0 + exp(-(Vₘ + 35.0) / 10.0)), name = :h)]

kdr_kinetics = [
    Gate(AlphaBeta,
         ifelse(Vₘ == -55.0, 0.1, (0.01 * (Vₘ + 55.0)) / (1.0 - exp(-(Vₘ + 55.0) / 10.0))),
         0.125 * exp(-(Vₘ + 65.0) / 80.0),
         p = 4, name = :n)]

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS / cm^2)
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS / cm^2)
@named leak = IonChannel(Leak, max_g = 0.3mS / cm^2)

channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV])

@named holding_current = Bias(5000pA);

dynamics_1 = HodgkinHuxley(channels, reversals);
dynamics_2 = HodgkinHuxley(channels, reversals);

@named neuron1 = Compartment(Vₘ, dynamics_1;
                             geometry = Cylinder(radius = 25µm, height = 400µm),
                             stimuli = [holding_current]);
@named neuron2 = Compartment(Vₘ, dynamics_2;
                            geometry = Cylinder(radius = 25µm, height = 400µm))

# Event-based Synaptic model
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@variables m(t) = 0.0
@parameters τsyn = 20 # 20 ms
syn_kinetics = Gate(SimpleGate, m, [D(m) ~ -m/τsyn])
event_model = ConstantValueEvent(1.0, m) # increment `m` gate by 1 for each incoming AP
@named Glut = SynapticChannel(event_model, Cation, [syn_kinetics]; max_s = 30nS);

topology = NetworkTopology([neuron1, neuron2], [Glut]);

topology[neuron1, neuron2] = Glut
reversal_map = Dict([Glut => EGlut])

@named net = NeuronalNetworkSystem(topology, reversal_map)
total_time = 250.0
sim = Simulation(net, time = total_time * ms)
solution = solve(sim, Rosenbrock23(), abstol = 1e-3, reltol = 1e-3, saveat = 0.2)

# Plot at 5kHz sampling
plot(solution(0.0:0.2:total_time, idxs = [neuron1.Vₘ, neuron2.Vₘ]))
