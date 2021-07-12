cd("Conductor.jl/test/")
include(joinpath(@__DIR__, "prinz_kinetics.jl"))
include(joinpath(@__DIR__, "prinz_synapses.jl"))

using OrdinaryDiffEq, Plots

# Figure 3e Prinz (rough attempt)
# AB/PD 3 (Note: I think this is supposed to be ABPD 2, but that one results in silent network
@named ABPD = Soma([NaV( 200mS/cm^2),
                    CaT( 2.5mS/cm^2),
                    CaS(   4mS/cm^2),
                    KA(   50mS/cm^2),
                    KCa(   5mS/cm^2),
                    Kdr( 50mS/cm^2),
                    H(   .01mS/cm^2),
                    leak(  0mS/cm^2)],
                    gradients, area = area, V0 = -50mV, aux = [calcium_conversion]);

# LP 4
@named LP = Soma([NaV( 100mS/cm^2),
                  CaT(   0mS/cm^2),
                  CaS(   4mS/cm^2),
                  KA(   20mS/cm^2),
                  KCa(   0mS/cm^2),
                  Kdr(  25mS/cm^2),
                  H(   .05mS/cm^2),
                  leak(.03mS/cm^2)],
                  gradients, area = area, V0 = -50mV, aux = [calcium_conversion]);

# PY 1
@named PY = Soma([NaV( 100mS/cm^2),
                  CaT( 2.5mS/cm^2),
                  CaS(   2mS/cm^2),
                  KA(   50mS/cm^2),
                  KCa(   0mS/cm^2),
                  Kdr( 125mS/cm^2),
                  H(   .05mS/cm^2),
                  leak(.01mS/cm^2)],
                  gradients, area = area, V0 = -50mV, aux = [calcium_conversion]);

topology = [ABPD => (LP, Glut(30nS)),
            ABPD => (LP, Chol(30nS)),
            ABPD => (PY, Glut(10nS)),
            ABPD => (PY, Chol(3nS)),
            LP   => (ABPD, Glut(30nS)),
            LP   => (PY, Glut(1nS)),
            PY   => (LP, Glut(30nS))];

network = Network([ABPD, LP, PY], topology)

t = 7500
sim = Simulation(network, time = t*ms)
solution = solve(sim, Rosenbrock23())

# Plot at 5kHz sampling
plot(solution; plotdensity=Int(t*5), vars = [ABPD.sys.Vₘ, LP.sys.Vₘ, PY.sys.Vₘ])

