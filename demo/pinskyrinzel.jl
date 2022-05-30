
cd("/home/wsphil/git/Conductor.jl/demo")
include("traub_kinetics.jl")

using ModelingToolkit, Unitful

import Unitful: µF, cm, pA, µA, nA

@parameters ϕ = 0.013 β = 0.075
@named calcium_conversion = ODESystem([D(Caᵢ) ~ -ϕ*ICa - β*Caᵢ]);


reversals = Equilibria(Pair[Sodium    =>  120.0mV,
                            Potassium =>  -15.0mV,
                            Leak      =>    0.0mV,
                            Calcium   =>  140.0mV]);

capacitance = 3µF/cm^2

# Pinsky modifies NaV to have instantaneous activation,
# so we can ignore tau
pinsky_nav_kinetics = [convert(Gate{SteadyState}, nav_kinetics[1]),
                       nav_kinetics[2]]

@named NaV = IonChannel(Sodium, pinsky_nav_kinetics) 

@named Iₑ = IonCurrent(NonIonic)
electrode_pulse = Iₑ ~ IfElse.ifelse(t > 100.0, IfElse.ifelse(t < 200.0, ustrip(Float64, µA, 5000nA), 0.0), 0.0)

@named soma = Compartment(Vₘ, [NaV(30mS/cm^2), Kdr(15mS/cm^2), leak(0.1mS/cm^2)], reversals[1:3], capacitance = capacitance)

@named dendrite = Compartment(Vₘ, [KAHP(0.8mS/cm^2), CaS(10mS/cm^2), KCa(15mS/cm^2), leak(0.1mS/cm^2)], reversals[2:4], capacitance = capacitance, extensions = [calcium_conversion])

@named g_c = AxialConductance([Gate(ConstantValue, 0.5, name = :p)], max_g = 2.1mS/cm^2)
jxn = Junction(soma => dendrite, g_c, symmetric = true)

test = MultiCompartment([jxn])


sim = Simulation(test, time=5000ms)

sol = solve(sim, Rosenbrock23())

plot(sol)
