αₘ(V) = ifelse(V == -40.0, 1.0, (0.1*(V + 40.0))/(1.0 - exp(-(V + 40.0)/10.0)))
βₘ(V) = 4.0*exp(-(V + 65.0)/18.0)
αₕ(V) = 0.07*exp(-(V + 65.0)/20.0)
βₕ(V) = 1.0/(1.0 + exp(-(V + 35.0)/10.0))
αₙ(V) = ifelse(V == -55.0, 0.1, (0.01*(V + 55.0))/(1.0 - exp(-(V + 55.0)/10.0)))
βₙ(V) = 0.125 * exp(-(V + 65.0)/80.0)

@register_symbolic αₘ(V)
@register_symbolic βₘ(V)
@register_symbolic αₕ(V)
@register_symbolic βₕ(V)
@register_symbolic αₙ(V)
@register_symbolic βₙ(V)

ModelingToolkit.get_unit(op::typeof(αₘ), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(βₘ), args) = NoUnits
ModelingToolkit.get_unit(αₕ::typeof(αₕ), args) = NoUnits
ModelingToolkit.get_unit(βₕ::typeof(βₕ), args) = NoUnits
ModelingToolkit.get_unit(αₙ::typeof(αₙ), args) = NoUnits
ModelingToolkit.get_unit(βₙ::typeof(βₙ), args) = NoUnits

