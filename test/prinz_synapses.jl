# Synaptic kinetics
syn∞ = 1/(1 + exp((-35 - Vₘ)/5))
glut_kinetics = Gate(SteadyStateTau, :s,
                     syn∞,
                     (1 - syn∞)/(1/40), 1)

chol_kinetics = Gate(SteadyStateTau, :s,
                     syn∞,
                     (1 - syn∞)/(1/100), 1)

EGlut = Equilibrium{Leak}(-70mV, :Glut) # NOTE: -70mV reversal => IPSP
EChol = Equilibrium{Leak}(-80mV, :Chol) # Leak as alias for non-specific ion current

@named Glut = SynapticChannel(Leak, [glut_kinetics], EGlut);
@named Chol = SynapticChannel(Leak, [chol_kinetics], EChol);


