
using LinearAlgebra, OrdinaryDiffEq, ModelingToolkit, Graphs
LinearAlgebra.BLAS.set_num_threads(6)

include(joinpath(@__DIR__, "traub_kinetics.jl"))
include(joinpath(@__DIR__, "pinsky_setup.jl"))

# Base dynamics
@named I_s_holding = Bias(-0.5µA / 0.5)
@named I_d_holding = Bias(0.0µA / (1 - 0.5))

soma_dynamics = HodgkinHuxley([NaV(30mS / cm^2),
                               Kdr(15mS / cm^2),
                               leak(0.1mS / cm^2)],
                               reversals[1:3]);

@named soma = Compartment(Vₘ, soma_dynamics;
                          geometry = Unitless(0.5), # FIXME: should take p param
                          capacitance = capacitance,
                          stimuli = [I_s_holding])

dendrite_dynamics = HodgkinHuxley([KAHP(0.8mS / cm^2),
                                   CaS(10mS / cm^2),
                                   KCa(15mS / cm^2),
                                   leak(0.1mS / cm^2)],
                                  reversals[2:4]);

@named dendrite = Compartment(Vₘ, dendrite_dynamics;
                              extensions = [calcium_conversion],
                              geometry = Unitless(0.5), # FIXME: should take p param
                              capacitance = capacitance,
                              stimuli = [I_d_holding])

@named gc_soma = AxialConductance([Gate(inv(p), name = :ps)], max_g = gc_val)
@named gc_dendrite = AxialConductance([Gate(inv(1 - p), name = :pd)], max_g = gc_val)
topology = MultiCompartmentTopology([soma, dendrite]);
add_junction!(topology, soma, dendrite, (gc_soma, gc_dendrite))
@named mcneuron = MultiCompartment(topology)

###########################################################################################
# Synapse models
###########################################################################################
import Conductor: NMDA, AMPA

@parameters τNMDA = 150.0 τAMPA = 2.0
@variables S(t) = 0.0
@named NMDAChan = SynapticChannel(ConstantValueEvent(S; saturation = 125.0), NMDA,
                                  [Gate(inv(1 + 0.28 * exp(-0.062(Vₘ - 60.0))); name = :e),
                                   Gate(S, [D(S) ~ -S/τNMDA]),
                                   Gate(inv(1 - p), name = :pnmda)],
                                  max_s = 0.02mS)

@named AMPAChan = SynapticChannel(ConstantValueEvent(S; threshold = 20mV), AMPA,
                                  [Gate(S, [D(S) ~ -S/τAMPA]),
                                   Gate(inv(1 - p), name = :pampa)],
                                  max_s = 0.04mS)

ENMDA = EquilibriumPotential(NMDA, 60mV, name = :NMDA)
EAMPA = EquilibriumPotential(AMPA, 60mV, name = :AMPA)
revmap = Dict([NMDAChan => ENMDA, AMPAChan => EAMPA])

# Create a neuron group
n_neurons = 100
@named neurons = Population(mcneuron, n_neurons);

# A brief stimulation to trigger the network
@named I_pulse = PulseTrain(amplitude = 5.0µA / 0.5,
                            duration = 50ms,
                            delay = 150.0ms,
                            offset = -0.5µA / 0.5)

# apply the stimulus to the root compartment of the 4th neuron
add_stimuli!(neurons, [I_pulse], (4,1))
topology = NetworkTopology(neurons, [NMDAChan, AMPAChan]);

# Generate connectivity graphs
nmda_g = random_regular_digraph(n_neurons, fld(n_neurons, 5), dir = :in)
ampa_g = random_regular_digraph(n_neurons, fld(n_neurons, 5), dir = :in)

# Use graphs to create synapses from soma (output) to dendrites (inputs) for neurons
for (i, e) in enumerate(edges(nmda_g))
    add_synapse!(topology, src(e), (dst(e),2),
                 NMDAChan, 1.0)
end

for (i, e) in enumerate(edges(ampa_g))
    add_synapse!(topology, src(e), (dst(e),2),
                 AMPAChan, 1.0)
end

@named net = NeuronalNetworkSystem(topology, revmap);
t_total = 2000.0
prob = Simulation(net, t_total*ms)
@time sol = solve(prob, RK4(); adaptive = false, dt = 0.05);

############################################################################################
# Pinsky and Rinzel displayed their results as a plot of N neurons over 20mV
using Statistics, Plots

# helper function
indexof(sym, syms) = findfirst(isequal(sym), syms)

simp = ODESystem(net)
dvs = states(simp)
interpolated = sol(1:0.2:t_total,
                   idxs = [indexof(x.soma.Vₘ, dvs) for x in neurons[5:end]])
abovethold = reduce(hcat, interpolated.u) .> 10.0
final = mean(abovethold, dims = 1)'
plot(final)

