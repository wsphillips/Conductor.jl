
include("traub_kinetics.jl")

ModelingToolkit.setdefault(ϕ, 0.13)
gradients = Equilibria(Pair[Sodium    =>  120.0mV,
                            Potassium =>  -15.0mV,
                            Leak      =>    0.0mV,
                            Calcium   =>  140.0mV]);

capacitance = 3µF/cm^2

# convert NaV to have instantaneous activation...



