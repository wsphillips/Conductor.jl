import Unitful: µF, pA, µA, nA, µS
@parameters ϕ = 0.13 β = 0.075 p = 0.5
@named calcium_conversion = ODESystem([D(Caᵢ) ~ -ϕ*ICa - β*Caᵢ];
                                      checks = ModelingToolkit.CheckComponents);
reversals = Equilibria([Sodium    =>  120.0mV,
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

