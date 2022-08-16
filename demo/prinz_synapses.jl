import Conductor: Cholinergic, Glutamatergic
# Synaptic kinetics
Vₓ = ExtrinsicPotential()
syn∞(V) = 1/(1 + exp((-35 - V)/5))
Glut₊τsyn(V) = (1 - syn∞(V))/(1/40)
Chol₊τsyn(V) = (1 - syn∞(V))/(1/100)

@register_symbolic syn∞(V)
@register_symbolic Glut₊τsyn(V)
@register_symbolic Chol₊τsyn(V)

ModelingToolkit.get_unit(op::typeof(syn∞), args) = NoUnits
ModelingToolkit.get_unit(op::typeof(Glut₊τsyn), args) = ms
ModelingToolkit.get_unit(op::typeof(Chol₊τsyn), args) = ms

glut_kinetics = Gate(SteadyStateTau, syn∞(Vₓ), Glut₊τsyn(Vₓ), name = :z1)
chol_kinetics = Gate(SteadyStateTau, syn∞(Vₓ), Chol₊τsyn(Vₓ), name = :z2)

EGlut = Equilibrium(Glutamatergic, -70mV, name = :Glut) # NOTE: -70mV reversal => IPSP
EChol = Equilibrium(Cholinergic, -80mV, name = :Chol) # Leak as alias for non-specific ion current

@named Glut = SynapticChannel(Glutamatergic, [glut_kinetics]);
@named Chol = SynapticChannel(Cholinergic, [chol_kinetics]);


