
include(joinpath(@__DIR__, "traub_kinetics.jl")

import Unitful: µF, pA, µA, nA

@parameters ϕ = 0.13 β = 0.075
@named calcium_conversion = ODESystem([D(Caᵢ) ~ -ϕ*ICa - β*Caᵢ]);

reversals = Equilibria(Pair[Sodium    =>  120.0mV,
                            Potassium =>  -15.0mV,
                            Leak      =>    0.0mV,
                            Calcium   =>  140.0mV]);

capacitance = 3.0µF/cm^2

# Pinsky modifies NaV to have instantaneous activation, so we can ignore tau
pinsky_nav_kinetics = [convert(Gate{SteadyState}, nav_kinetics[1]), nav_kinetics[2]]
@named NaV = IonChannel(Sodium, pinsky_nav_kinetics) 

# No inactivation term for calcium current in Pinsky model
pinsky_ca_kinetics = [ca_kinetics[1]]
@named CaS = IonChannel(Calcium, pinsky_ca_kinetics)

@named Iₛ = IonCurrent(NonIonic)
soma_holding = Iₛ ~ ustrip(Float64, µA, 1.5µA)

@named I_d = IonCurrent(NonIonic)
dendrite_holding = I_d ~ ustrip(Float64, µA, 0.0µA)

@named soma = Compartment(Vₘ,
                         [NaV(30mS/cm^2),
                          Kdr(15mS/cm^2),
                          leak(0.1mS/cm^2)],
                          reversals[1:3],
                          capacitance = capacitance,
                          stimuli = [soma_holding])

@named dendrite = Compartment(Vₘ,
                             [KAHP(0.8mS/cm^2),
                              CaS(10mS/cm^2),
                              KCa(15mS/cm^2),
                              leak(0.1mS/cm^2)],
                             reversals[2:4],
                             capacitance = capacitance,
                             extensions = [calcium_conversion],
                             stimuli = [dendrite_holding])

@named g_c = AxialConductance([Gate(ConstantValue, 2, name = :p)],
                              max_g = 2.1mS/cm^2)
jxn = Junction(soma => dendrite, g_c, symmetric = true)
mc_neuron = MultiCompartment([jxn])

simp = Simulation(mc_neuron, time=2000ms, return_system = true)
# Explicitly setting u0 for now...
prob = ODAEProblem(simp, [-4.6, 0.999, 0.001, 0.2, -4.5, 0.01, 0.009, .007], (0., 2000), [])

using OrdinaryDiffEq, Plots
# Pinsky & Rinzel originally solved using RK4 and dt=0.05
#sol = solve(prob, RK4(), dt=0.05, maxiters=1e9)
sol = solve(prob, Rosenbrock23(), reltol=1e-8, abstol=1e-8)
plot(sol, vars=[soma.Vₘ])
