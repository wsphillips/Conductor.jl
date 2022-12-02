# Example of synapses
using Conductor, OrdinaryDiffEq, Unitful, ModelingToolkit, Plots
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms, nS, pS
import Conductor: Na, K

############################################################################################
# Setup
############################################################################################

Vₘ = ParentScope(MembranePotential(-65mV));

nav_kinetics = [Gate(AlphaBeta,
                     ifelse(Vₘ == -40.0, 1.0,
                            (0.1 * (Vₘ + 40.0)) / (1.0 - exp(-(Vₘ + 40.0) / 10.0))),
                     4.0 * exp(-(Vₘ + 65.0) / 18.0),
                     p = 3, name = :m)
                Gate(AlphaBeta,
                     0.07 * exp(-(Vₘ + 65.0) / 20.0),
                     1.0 / (1.0 + exp(-(Vₘ + 35.0) / 10.0)), name = :h)];

kdr_kinetics = [Gate(AlphaBeta,
                     ifelse(Vₘ == -55.0, 0.1,
                            (0.01 * (Vₘ + 55.0)) / (1.0 - exp(-(Vₘ + 55.0) / 10.0))),
                     0.125 * exp(-(Vₘ + 65.0) / 80.0),
                     p = 4, name = :n)];

@named NaV = IonChannel(Sodium, nav_kinetics, max_g = 120mS / cm^2);
@named Kdr = IonChannel(Potassium, kdr_kinetics, max_g = 36mS / cm^2);
@named leak = IonChannel(Leak, max_g = 0.3mS / cm^2);

channels = [NaV, Kdr, leak];
reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV]);
dynamics = HodgkinHuxley(channels, reversals);

# Stimulate first neuron with a holding current to drive APs
@named holding_current = Bias(5000pA);
@named neuron1 = Compartment(Vₘ, dynamics;
                             geometry = Cylinder(radius = 25µm, height = 400µm),
                             stimuli = [holding_current]);

@named neuron2 = Compartment(Vₘ, dynamics;
                             geometry = Cylinder(radius = 25µm, height = 400µm));

############################################################################################
# Event-based Synaptic model
############################################################################################

EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@variables m(t) = 0.0
@parameters τsyn = 25 # 25 ms
syn_kinetics = Gate(SimpleGate, m, [D(m) ~ -m/τsyn])
event_model = ConstantValueEvent(1.0, m) # increment `m` gate by 1 for each incoming AP
@named ExpAMPA = SynapticChannel(event_model, Cation, [syn_kinetics]; max_s = 30nS)

top = NetworkTopology([neuron1, neuron2], [ExpAMPA]);

top[neuron1, neuron2] = ExpAMPA(10nS)
rev_map = Dict([ExpAMPA => EGlut])

@named net = NeuronalNetworkSystem(top, rev_map)

total_time = 250.0
sim = Simulation(net, total_time * ms)
sol = solve(sim, Rosenbrock23(), abstol=1e-3, reltol=1e-3);

plot(plot(sol, idxs = [neuron1.Vₘ]),
     plot(sol, idxs = [neuron2.Vₘ]),
     layout=(2,1))

############################################################################################
# Continuously integrated synaptic model
############################################################################################
EGlut = Equilibrium(Cation, 0mV, name = :Glut)
Vₓ = ExtrinsicPotential()
syn∞ = 1 / (1 + exp((-35 - Vₓ) / 5))
tausyn = (1 - syn∞) / (1 / 40)
syn_kinetics2 = Gate(SteadyStateTau, syn∞, tausyn, name = :z)

@named IntAMPA = SynapticChannel(IntegratedSynapse(), Cation, [syn_kinetics2]; max_s = 30nS)

top2 = NetworkTopology([neuron1, neuron2], [IntAMPA]);

top2[neuron1, neuron2] = IntAMPA(30nS)
rev_map2 = Dict([IntAMPA => EGlut])

@named net2 = NeuronalNetworkSystem(top2, rev_map2)

total_time = 250.0
sim2 = Simulation(net2, total_time * ms)
sol2 = solve(sim2, Rosenbrock23(), abstol = 1e-3, reltol = 1e-3);

plot(plot(sol2, idxs = [neuron1.Vₘ]),
     plot(sol2, idxs = [neuron2.Vₘ]),
     layout=(2,1))



