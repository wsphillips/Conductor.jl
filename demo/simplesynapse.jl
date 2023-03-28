# Example of synapses
using Conductor, Unitful, ModelingToolkit, OrdinaryDiffEq, Plots
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

reversals = Equilibria([Na => 50.0mV, K => -77.0mV, Leak => -54.4mV]);
dynamics = HodgkinHuxley([NaV, Kdr, leak], reversals);
geometry = Cylinder(radius = 25µm, height = 400µm)

# Stimulate first neuron with a holding current to drive APs
@named holding_current = Bias(5000pA);
@named neurons = Population(Compartment(Vₘ, dynamics; geometry), 2);

Conductor.add_stimuli!(neurons, [holding_current], 1)

############################################################################################
# Event-based Synaptic model
############################################################################################

EGlut = Equilibrium(Cation, 0mV, name = :Glut)
@variables m(t) = 0.0
@parameters τsyn = 25 # 25 ms
syn_kinetics = Gate(SimpleGate, m, [D(m) ~ -m/τsyn])
event_model = ConstantValueEvent(m) # increment `m` gate by 1 for each incoming AP
@named ExpAMPA = SynapticChannel(event_model, Cation, [syn_kinetics]; max_s = 30nS)

top = NetworkTopology(neurons, [ExpAMPA]);

top[neurons[1], neurons[2]] = ExpAMPA(10nS)
rev_map = Dict([ExpAMPA => EGlut])

@named net = NeuronalNetworkSystem(top, rev_map)

total_time = 250.0
sim = Simulation(net, total_time * ms)
sol = solve(sim, Rosenbrock23(), abstol=1e-3, reltol=1e-3);

plot(plot(sol, idxs = [neurons[1].Vₘ]),
     plot(sol, idxs = [neurons[2].Vₘ]),
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

top2 = NetworkTopology(neurons, [IntAMPA]);

top2[neurons[1], neurons[2]] = IntAMPA(300nS)
rev_map2 = Dict([IntAMPA => EGlut])

@named net2 = NeuronalNetworkSystem(top2, rev_map2)

total_time = 250.0
sim2 = Simulation(net2, (0.0ms, total_time * ms))
sol2 = solve(sim2, Rosenbrock23(), abstol = 1e-3, reltol = 1e-3);

plot(plot(sol2, idxs = [neurons[1].Vₘ]),
     plot(sol2, idxs = [neurons[2].Vₘ]),
     layout=(2,1))

