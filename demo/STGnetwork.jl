# Synaptically connected STG model network built as in Prinz et al Nature Neurosci 2004

using OrdinaryDiffEq, Plots
import Unitful.nS
include(joinpath(@__DIR__, "prinz_kinetics.jl"));
include(joinpath(@__DIR__, "prinz_synapses.jl"));

# Figure 3e Prinz (qualitatively close to expected output)
# AB/PD 2 
@named ABPD = Compartment(Vₘ, [NaV(100mS/cm^2),
                    CaT(2.5mS/cm^2),
                    CaS(6mS/cm^2),
                    KA(50mS/cm^2),
                    KCa(5mS/cm^2),
                    Kdr(100mS/cm^2),
                    H(.01mS/cm^2),
                    leak(0mS/cm^2)],
                    gradients, geometry = geo, extensions = [calcium_conversion]);

ABPD # display

# LP 4
@named LP = Compartment(Vₘ, [NaV( 100mS/cm^2),
                  CaT(   0mS/cm^2),
                  CaS(   4mS/cm^2),
                  KA(   20mS/cm^2),
                  KCa(   0mS/cm^2),
                  Kdr(  25mS/cm^2),
                  H(   .05mS/cm^2),
                  leak(.03mS/cm^2)],
                  gradients, geometry = geo, extensions = [calcium_conversion]);

LP # display

# PY 1
@named PY = Compartment(Vₘ, [NaV( 100mS/cm^2),
                  CaT( 2.5mS/cm^2),
                  CaS(   2mS/cm^2),
                  KA(   50mS/cm^2),
                  KCa(   0mS/cm^2),
                  Kdr( 125mS/cm^2),
                  H(   .05mS/cm^2),
                  leak(.01mS/cm^2)],
                  gradients, geometry = geo, extensions = [calcium_conversion]);

PY # display

topology = [Synapse(ABPD => LP, Glut(30nS), EGlut),
            Synapse(ABPD => LP, Chol(30nS), EChol),
            Synapse(ABPD => PY, Glut(10nS), EGlut),
            Synapse(ABPD => PY, Chol(3nS), EChol),
            Synapse(LP   => ABPD, Glut(30nS), EGlut),
            Synapse(LP   => PY, Glut(1nS), EGlut),
            Synapse(PY   => LP, Glut(30nS), EGlut)];

network = NeuronalNetworkSystem(topology)

sim = Simulation(network, time = 10000ms)
solution = solve(sim, Rosenbrock23())

# Plot at 5kHz sampling
fig = plot(solution; plotdensity=Int(10000*5), size=(1200,800), vars = [ABPD.Vₘ, LP.Vₘ, PY.Vₘ])
fig

# Uncomment and eval `png(...)` to save as PNG
# png(fig, "figure3e_simulated")

