# Gates

Gates are the most basic building block for describing dynamics in Conductor.jl. A gate
integrates zero or more inputs into a single unitless weighting value. If desired, we can
explicitly define the symbolic expression that models the output of a gate:

```julia
# Voltage-dependent sigmoid activation
Vₘ = MembranePotential()
@named sigmoid = Gate(SimpleGate, inv(1 + exp(-Vₘ)))
```

The dynamics of a `Gate{SimpleGate}` are just an algebraic equation:

```julia
Conductor.get_eqs(sigmoid, nothing) # internal API call
    # 1-element Vector{Symbolics.Equation}:
    # sigmoid(t) ~ 1 / (1 + exp(-Vₘ(t)))
```

Alternatively, we can take advantage of customized constructors, which let us describe
gate dynamics in a more canonical, domain-friendly form. For instance, in the classic
Hodgkin-Huxley formalism, gates are analogous to the "gating particles" used to describe the
kinetics of voltage-gated ion channels. Ionic channel kinetics are often described in terms
of forward (``\alpha``) and reverse (``\beta``) reaction rates:

```julia

# Sodium channel inactivation kinetics, defined by forward and reverse rxn rates
@named h = Gate(
    AlphaBeta,
    # the α (forward) reaction rate
    0.07*exp(-(Vₘ+65.0)/20.0),
    # the β (reverse) reaction rate
    1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0))
)
```

Now when we ask for the dynamics of a `Gate{AlphaBeta}`, we get the differential form:
```julia
Conductor.get_eqs(h, nothing)
#1-element Vector{Symbolics.Equation}:
#Differential(t)(h(t)) ~ ((0.07exp(-3.25 - 0.05Vₘ(t))) / ( 1.0 / (1.0 + exp(-3.5 - 0.1Vₘ(t)))
#+ 0.07exp(-3.25 - 0.05V ₘ(t))) - h(t))*(1.0 / (1.0 + exp(-3.5 - 0.1Vₘ(t))) + 0.07e xp(-3.25
#- 0.05Vₘ(t)))
```
...which is the equivalent to:
```math
\frac{dh}{dt} = \alpha_h (1-n)-\beta_h n
```

Customized gates also enable dynamic equation writing. In some cases, the dynamics of a gate
may depend on the structure or properties of the parent `CompartmentSystem`.

```math
S\prime={\displaystyle \sum_{j}H(V_{s,j}-10)-S/150} 
```

Multiple gates can be combined as a product that represents the state-dependent activation
of a conductance.

