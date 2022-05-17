import Conductor: Cholinergic, Glutamatergic
# Synaptic kinetics
syn∞ = 1/(1 + exp((-35 - Vₘ)/5))
glut_kinetics = Gate(SteadyStateTau,
                     syn∞,
                     (1 - syn∞)/(1/40), name = :m)

chol_kinetics = Gate(SteadyStateTau,
                     syn∞,
                     (1 - syn∞)/(1/100), name = :m)

EGlut = Equilibrium(Glutamatergic, -70mV, name = :Glut) # NOTE: -70mV reversal => IPSP
EChol = Equilibrium(Cholinergic, -80mV, name = :Chol) # Leak as alias for non-specific ion current

@named Glut = SynapticChannel(Glutamatergic, [glut_kinetics]);
@named Chol = SynapticChannel(Cholinergic, [chol_kinetics]);


