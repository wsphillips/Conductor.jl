
include(joinpath(@__DIR__, "traub_kinetics.jl"))

using LinearAlgebra, OrdinaryDiffEq
LinearAlgebra.BLAS.set_num_threads(6)
import Unitful: µF, pA, µA, nA, µS

@parameters ϕ = 0.13 β = 0.075 p = 0.5

@named calcium_conversion = ODESystem([D(Caᵢ) ~ -ϕ*ICa - β*Caᵢ]);

reversals = Equilibria(Pair[Sodium    =>  120.0mV,
                            Potassium =>  -15.0mV,
                            Leak      =>    0.0mV,
                            Calcium   =>  140.0mV]);

capacitance = 3.0µF/cm^2
gc_val = 2.1mS/cm^2

# Pinsky modifies NaV to have instantaneous activation, so we can ignore tau
pinsky_nav_kinetics = [convert(Gate{SimpleGate}, nav_kinetics[1]), nav_kinetics[2]]
@named NaV = IonChannel(Sodium, pinsky_nav_kinetics) 

# No inactivation term for calcium current in Pinsky model
pinsky_ca_kinetics = [ca_kinetics[1]]
@named CaS = IonChannel(Calcium, pinsky_ca_kinetics)

is_val = ustrip(Float64, µA, -0.20µA)/p
@named Iₛ = IonCurrent(NonIonic, is_val, dynamic = false)
soma_holding = Iₛ ~ is_val

id_val = ustrip(Float64, µA, 0.0µA)/(1-p)
@named I_d = IonCurrent(NonIonic, id_val, dynamic = false)
dendrite_holding = I_d ~ id_val

soma_dynamics = HodgkinHuxley(Vₘ,
                         [NaV(30mS/cm^2),
                          Kdr(15mS/cm^2),
                          leak(0.1mS/cm^2)],
                          reversals[1:3];
                          geometry = Unitless(0.5), # FIXME: should take p param
                          capacitance = capacitance,
                          stimuli = [soma_holding])

@named soma = Compartment(soma_dynamics)

dendrite_dynamics = HodgkinHuxley(Vₘ,
                             [KAHP(0.8mS/cm^2),
                              CaS(10mS/cm^2),
                              KCa(15mS/cm^2),
                              leak(0.1mS/cm^2)],
                             reversals[2:4],
                             geometry = Unitless(0.5), # FIXME: should take p param
                             capacitance = capacitance,
                             stimuli = [dendrite_holding])

@named dendrite = Compartment(dendrite_dynamics,
                              extensions = [calcium_conversion])


@named gc_soma = AxialConductance([Gate(SimpleGate, inv(p), name = :ps)],
                                  max_g = gc_val)
@named gc_dendrite = AxialConductance([Gate(SimpleGate, inv(1-p), name = :pd)],
                                      max_g = gc_val)

topology = Conductor.MultiCompartmentTopology([soma, dendrite]);

Conductor.add_junction!(topology, soma,  dendrite, gc_soma, symmetric = false)
Conductor.add_junction!(topology, dendrite,  soma, gc_dendrite, symmetric = false)

@named mcneuron = MultiCompartment(topology)

###########################################################################################
# Synapse models
###########################################################################################
import Conductor: NMDA, AMPA, HeavisideSum

@named NMDAChan = SynapticChannel(NMDA,
                [Gate(SimpleGate, inv(1 + 0.28*exp(-0.062(Vₘ - 60.))); name = :e),
                 Gate(HeavisideSum; threshold = 10mV, decay = 150, saturation = 125, name = :S),
                 Gate(SimpleGate, inv(1-p), name = :pnmda)],
                 max_s = 0.014mS, aggregate = true)

@named AMPAChan = SynapticChannel(AMPA,
                                [Gate(HeavisideSum, threshold = 20mV, decay = 2;
                                      name = :u),
                                 Gate(SimpleGate, inv(1-p), name = :pampa)],
                                 max_s = 0.0045mS, aggregate = true)

ENMDA = EquilibriumPotential(NMDA, 60mV, name = :NMDA)
EAMPA = EquilibriumPotential(AMPA, 60mV, name = :AMPA)
revmap = Dict([NMDAChan => ENMDA, AMPAChan => EAMPA])

# A single neuron was "briefly" stimulated to trigger the network
#is_stim_val = IfElse.ifelse(t < 50.0, ustrip(Float64, µA, 2.0µA)/p, ustrip(Float64, µA, -0.2µA)/p)
is_stim_val = ustrip(Float64, µA, 1.0µA)/p
@named Istim = IonCurrent(NonIonic, is_stim_val)
soma_stim = Istim ~ is_stim_val

soma_stim_dynamics = HodgkinHuxley(Vₘ,
                                   [NaV(30mS/cm^2),
                                    Kdr(15mS/cm^2),
                                    leak(0.1mS/cm^2)],
                                    reversals[1:3];
                                    geometry = Unitless(0.5), # FIXME: should take p param
                                    capacitance = capacitance,
                                    stimuli = [soma_stim])

@named soma_stimulated = Compartment(soma_stim_dynamics)

mcstim_topology = Conductor.MultiCompartmentTopology([soma_stimulated, dendrite]);

Conductor.add_junction!(mcstim_topology, soma_stimulated,  dendrite, gc_soma, symmetric = false)
Conductor.add_junction!(mcstim_topology, dendrite,  soma_stimulated, gc_dendrite, symmetric = false)

@named mcneuron_stim = MultiCompartment(mcstim_topology)

# Need to introduce 10% gca variance as per Pinsky/Rinzel
neuronpopulation = [Conductor.replicate(mcneuron) for _ in 1:100];
neuronpopulation[4] = mcneuron_stim
topology = NetworkTopology(neuronpopulation, [NMDAChan, AMPAChan]);

using Graphs
nmda_g = random_regular_digraph(100, 20, dir=:in)
ampa_g = random_regular_digraph(100, 20, dir=:in)

for (i, e) in enumerate(edges(nmda_g))
    if src(e) == 4
        add_synapse!(topology, neuronpopulation[src(e)].soma_stimulated, neuronpopulation[dst(e)].dendrite, NMDAChan)
    else
        add_synapse!(topology, neuronpopulation[src(e)].soma, neuronpopulation[dst(e)].dendrite, NMDAChan)
    end
end
for (i, e) in enumerate(edges(ampa_g))
    if src(e) == 4
        add_synapse!(topology, neuronpopulation[src(e)].soma_stimulated, neuronpopulation[dst(e)].dendrite, AMPAChan)
    else
        add_synapse!(topology, neuronpopulation[src(e)].soma, neuronpopulation[dst(e)].dendrite, AMPAChan)
    end
end

@named net = NeuronalNetworkSystem(topology, revmap);
simp = Simulation(net, time = 2000ms, return_system = true)
prob = Simulation(net, time = 2000ms)

# on the first run, over 40% of the time is single threaded compilation time. Without
# compilation, solving takes roughly 200 seconds (6-8 cores) for 5 seconds of tspan.
@time sol = solve(prob, RadauIIA5(), abstol=1e-6, reltol=1e-6, saveat=0.2);
#plot(sol, vars = [x.soma.Vₘ for x in neuronpopulation])

# Pinsky and Rinzel displayed their results as a plot of N neurons over 20mV
indexof(sym,syms) = findfirst(isequal(sym),syms)
dvs = states(simp)
interpolated = sol(1:2000, idxs=[indexof(x.soma.Vₘ, dvs) for x in neuronpopulation[5:end]])
abovethold = reduce(hcat, interpolated.u) .> 20.0
using Statistics
final = mean(abovethold, dims=1)'

using Plots

plot(final)
