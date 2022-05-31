
cd("/home/wikphi/git/Conductor.jl/demo")
include("traub_kinetics.jl")

using OrdinaryDiffEq, Plots
import Unitful: µF, pA, µA, nA

@parameters ϕ = 0.13 β = 0.075
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
pinsky_ca_kinetics = [ca_kinetics[1]]
@named NaV = IonChannel(Sodium, pinsky_nav_kinetics) 
@named CaS = IonChannel(Calcium, pinsky_ca_kinetics)

@named Iₛ = IonCurrent(NonIonic)
soma_holding = Iₛ ~ ustrip(Float64, µA, 10µA)

@named I_d = IonCurrent(NonIonic)
dendrite_holding = I_d ~ ustrip(Float64, µA, 0.0µA)

@named soma = Compartment(Vₘ, [NaV(15mS/cm^2), Kdr(7.5mS/cm^2), leak(0.05mS/cm^2)], reversals[1:3], capacitance = capacitance, stimuli = [soma_holding])

@named dendrite = Compartment(Vₘ, [KAHP(0.4mS/cm^2), CaS(5mS/cm^2), KCa(7.5mS/cm^2), leak(0.05mS/cm^2)], reversals[2:4], capacitance = capacitance, extensions = [calcium_conversion], stimuli = [dendrite_holding])

@named g_c = AxialConductance([Gate(ConstantValue, 2, name = :p)], max_g = 2.1mS/cm^2)
jxn = Junction(soma => dendrite, g_c, symmetric = true)
#
test = MultiCompartment([jxn])
sim = Simulation(test, time=2000ms, return_system = true)

#sim = Simulation(dendrite, time=2000ms, return_system = true)
#prob = ODAEProblem(sim, [-4.6, 0.999, 0.001], (0., 500), [])
prob = ODAEProblem(sim, [0.2, -4.5, 0.01, 0.009, 0.007,-4.6, 0.999, 0.001], (0., 2000), [])
#prob = ODAEProblem(sim, [0.2, -4.5, 0.01, 0.009, 0.007], (0.,2000), [])
#sol = solve(sim, RK4(), dt=0.05)
sol = solve(prob, Rosenbrock23())
plot(sol)
