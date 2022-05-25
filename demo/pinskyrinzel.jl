
cd("/home/wsphil/git/Conductor.jl/demo")
include("traub_kinetics.jl")

using ModelingToolkit, Unitful

import Unitful: µF, cm

ModelingToolkit.setdefault(ϕ, 0.13)
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

@named soma = Compartment(Vₘ, [NaV, Kdr, leak], reversals[1:3], capacitance = capacitance)

@named dendrite = Compartment(Vₘ, [KAHP, CaS, KCa, leak], reversals[2:4], capacitance = capacitance)

Conductor.ScaledJunction()



