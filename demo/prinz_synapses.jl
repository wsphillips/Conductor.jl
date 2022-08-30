import Conductor: Cholinergic, Glutamatergic
# Synaptic kinetics
Vₓ = ExtrinsicPotential()
syn∞ = 1/(1 + exp((-35 - Vₓ)/5))
glut_kinetics = Gate(SteadyStateTau,
                     syn∞,
                     (1 - syn∞)/(1/40), name = :z1)

chol_kinetics = Gate(SteadyStateTau,
                     syn∞,
                     (1 - syn∞)/(1/100), name = :z2)

EGlut = Equilibrium(Glutamatergic, -70mV, name = :Glut) # NOTE: -70mV reversal => IPSP
EChol = Equilibrium(Cholinergic, -80mV, name = :Chol) # Leak as alias for non-specific ion current

@named Glut = SynapticChannel(Glutamatergic, [glut_kinetics]);
@named Chol = SynapticChannel(Cholinergic, [chol_kinetics]);

