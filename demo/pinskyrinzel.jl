using OrdinaryDiffEq, Plots

include(joinpath(@__DIR__, "traub_kinetics.jl"))
include(joinpath(@__DIR__, "pinsky_setup.jl"))

sim_time = 1500.0

# Compartment holding currents
@named Is_holding = Bias(-0.5µA / 0.5) # FIXME: we should be able to use 'p' parameter
@named Id_holding = Bias(0.0µA / (1 - 0.5))

soma_dynamics = HodgkinHuxley([NaV(30mS / cm^2), Kdr(15mS / cm^2), leak(0.1mS / cm^2)],
                              reversals[1:3]);

@named soma = Compartment(Vₘ, soma_dynamics;
                          capacitance = capacitance,
                          geometry = Unitless(0.5), # FIXME: support 'p' parameter value
                          stimuli = [Is_holding])

dendrite_dynamics = HodgkinHuxley([KAHP(0.8mS / cm^2), CaS(10mS / cm^2), KCa(15mS / cm^2),
                                   leak(0.1mS / cm^2)], reversals[2:4]);

@named dendrite = Compartment(Vₘ, dendrite_dynamics;
                              capacitance = capacitance,
                              geometry = Unitless(0.5),
                              stimuli = [Id_holding],
                              extensions = [calcium_conversion])

@named gc_soma = AxialConductance([Gate(inv(p), name = :ps)], max_g = gc_val)
@named gc_dendrite = AxialConductance([Gate(inv(1 - p), name = :pd)], max_g = gc_val)

topology = MultiCompartmentTopology([soma, dendrite]);
add_junction!(topology, soma, dendrite, (gc_soma, gc_dendrite))
@named mcneuron = MultiCompartment(topology);

# Note: Pinsky & Rinzel originally solved using RK4 and _fixed_ dt=0.05
# Here we let the solver use adaptive stepping, because its about 10X faster
# NOTE: To see activity, adjust the somatic/dendritic bias currents.
 prob = Simulation(mcneuron, sim_time*ms)
 sol = solve(prob, RK4())
 plot(sol, idxs=[soma.Vₘ], size=(1440,900))

# Published initial values
# prob = remake(prob; u0 = [-4.6, 0.999, 0.001, 0.2, -4.5, 0.01, 0.009, .007])

###########################################################################################
# Steady-state synaptic current
############################################################################################
import Conductor: NMDA
@named NMDAChan = IonChannel(NMDA,
                             [Gate(inv(1 + 0.28 * exp(-0.062(Vₘ - 60.0))); name = :e),
                              Gate(125.0*inv(1 - p), name = :pnmda)],
                              max_g = 4μS/cm^2)

ESyn = EquilibriumPotential(NMDA, 60mV, name = :syn)
stimulated_dendrite_dynamics = HodgkinHuxley([KAHP(0.8mS / cm^2),
                                              CaS(10mS / cm^2),
                                              KCa(15mS / cm^2),
                                              leak(0.1mS / cm^2),
                                              NMDAChan], [reversals[2:4]..., ESyn]);

@named stimulated_dendrite = Compartment(Vₘ, stimulated_dendrite_dynamics;
                              capacitance = capacitance,
                              geometry = Unitless(0.5),
                              stimuli = [Id_holding],
                              extensions = [calcium_conversion])

stimulated_topology = MultiCompartmentTopology([soma, stimulated_dendrite]);
add_junction!(stimulated_topology, soma, stimulated_dendrite, (gc_soma, gc_dendrite))
@named stimulated_mcneuron = MultiCompartment(stimulated_topology);

prob = Simulation(stimulated_mcneuron, sim_time * ms)
sol = solve(prob, RK4())
plot(sol, idxs = [stimulated_mcneuron.soma.Vₘ], size = (1440, 900))
