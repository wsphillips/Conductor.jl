# Gates

Gates are the most basic building block for describing dynamics in Conductor.jl. A gate
integrates 0+ inputs into a single unitless weighting value. If desired, we can explicitly
define the symbolic expression that models a gate:

```julia
# Voltage-dependent sigmoid activation
@named sigmoid = Gate(SimpleGate, inv(1 + exp(-Vₘ)))
```

Alternatively, we can take advantage of customized constructors, which let us describe
dynamics in a more canonical, domain-specific form. In the classic Hodgkin-Huxley formalism, gates are analogous to the "gating particles" used to describe the kinetics of voltage-gated ion channels. Their kinetics are often described in terms of forward (\alpha) and reverse (\beta) reaction rates:

```julia

Vₘ = MembranePotential()

# Sodium channel kinetics, described in terms for the forward and reverse rxn rates
nav_kinetics = [
    Gate(AlphaBeta,
         IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         4.0*exp(-(Vₘ + 65.0)/18.0), p = 3, name = :m)
    Gate(AlphaBeta,
         0.07*exp(-(Vₘ+65.0)/20.0),
         1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)), name = :h)]
```
Multiple gates can be combined
as a product that represents the state-dependent activation of a conductance.



