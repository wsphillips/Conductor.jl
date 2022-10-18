using OrdinaryDiffEq, Plots

include(joinpath(@__DIR__, "traub_kinetics.jl"))
include(joinpath(@__DIR__, "pinsky_setup.jl"))

sim_time = 1500.0

# Compartment holding currents
@named Is_holding = Bias(-0.5µA / 0.5) # FIXME: we should be able to use 'p' parameter
@named Id_holding = Bias(0.0µA / (1-0.5)) 

soma_dynamics = HodgkinHuxley(Vₘ,
                         [NaV(30mS/cm^2),
                          Kdr(15mS/cm^2),
                          leak(0.1mS/cm^2)],
                          reversals[1:3];
                          geometry = Unitless(0.5), # FIXME: support 'p' parameter value
                          capacitance = capacitance,
                          stimuli = [Is_holding]);

@named soma = Compartment(soma_dynamics)

dendrite_dynamics = HodgkinHuxley(Vₘ,
                             [KAHP(0.8mS/cm^2),
                              CaS(10mS/cm^2),
                              KCa(15mS/cm^2),
                              leak(0.1mS/cm^2)],
                             reversals[2:4],
                             geometry = Unitless(0.5),
                             capacitance = capacitance,
                             stimuli = [Id_holding]);

@named dendrite = Compartment(dendrite_dynamics, extensions = [calcium_conversion])

@named gc_soma = AxialConductance([Gate(SimpleGate, inv(p), name = :ps)], max_g = gc_val)
@named gc_dendrite = AxialConductance([Gate(SimpleGate, inv(1-p), name = :pd)], max_g = gc_val)

topology = MultiCompartmentTopology([soma,dendrite]);
add_junction!(topology, soma,  dendrite, (gc_soma, gc_dendrite))
@named mcneuron = MultiCompartment(topology);

# Note: Pinsky & Rinzel originally solved using RK4 and _fixed_ dt=0.05
# Here we let the solver use adaptive stepping, because its about 10X faster
# prob = Simulation(mcneuron, time=sim_time*ms)
# sol = solve(prob, RK4())
# plot(sol(0.0:0.2:sim_time, idxs=[soma.Vₘ]), size=(1440,900))

# Published initial values
# prob = remake(prob; u0 = [-4.6, 0.999, 0.001, 0.2, -4.5, 0.01, 0.009, .007])

###########################################################################################
# Steady synaptic inputs
############################################################################################
import Conductor: NMDA, AMPA, HeavisideSum

@named NMDAChan = SynapticChannel(NMDA,
                [Gate(SimpleGate, inv(1 + 0.28*exp(-0.062(Vₘ - 60.))); name = :e),
                 Gate(HeavisideSum; threshold = 10mV, decay = 150, saturation = 125, name = :S),
                 Gate(SimpleGate, inv(1-p), name = :pnmda)],
                 max_s = 0mS, aggregate = true)

# To simulate constant NMDA activation, we make a fake suprathreshold cell
artificial_Eleak = EquilibriumPotential(Leak, 20mV)
artificial_dynamics = HodgkinHuxley(Vₘ, [leak(1mS/cm^2)], [artificial_Eleak])
@named dummy = Compartment(artificial_dynamics)

ESyn = EquilibriumPotential(NMDA, 60mV, name = :syn)
topology = NetworkTopology([dummy, mcneuron], [NMDAChan(2µS)]);
topology[dummy, mcneuron.dendrite] = NMDAChan
revmap = [NMDAChan => ESyn]
network = NeuronalNetworkSystem(topology, revmap)

prob = Simulation(network, time=sim_time*ms)
sol = solve(prob, RK4())
plot(sol(0.0:0.2:sim_time, idxs=[mcneuron.soma.Vₘ]), size=(1440,900))

