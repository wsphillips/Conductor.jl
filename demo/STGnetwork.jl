# Synaptically connected STG model network built as in Prinz et al Nature Neurosci 2004

using OrdinaryDiffEq, Plots
import Unitful.nS
include(joinpath(@__DIR__, "prinz_kinetics.jl"));
include(joinpath(@__DIR__, "prinz_synapses.jl"));

# Figure 3e Prinz (qualitatively close to expected output)
# AB/PD 2 
ABPD_channels = [NaV(100mS/cm^2),
                 CaT(2.5mS/cm^2),
                 CaS(6mS/cm^2),
                 KA(50mS/cm^2),
                 KCa(5mS/cm^2),
                 Kdr(100mS/cm^2),
                 H(.01mS/cm^2),
                 leak(0mS/cm^2)];
ABPD_dynamics = HodgkinHuxley(Vₘ, ABPD_channels, gradients; geometry = geo);

@named ABPD = Compartment(ABPD_dynamics, extensions = [calcium_conversion]);

ABPD # display

# LP 4
LP_channels = [NaV( 100mS/cm^2),
               CaT(   0mS/cm^2),
               CaS(   4mS/cm^2),
               KA(   20mS/cm^2),
               KCa(   0mS/cm^2),
               Kdr(  25mS/cm^2),
               H(   .05mS/cm^2),
               leak(.03mS/cm^2)];
LP_dynamics = HodgkinHuxley(Vₘ, LP_channels, gradients; geometry = geo);
@named LP = Compartment(LP_dynamics, extensions = [calcium_conversion]);

LP # display

# PY 1
PY_channels = [NaV( 100mS/cm^2),
               CaT( 2.5mS/cm^2),
               CaS(   2mS/cm^2),
               KA(   50mS/cm^2),
               KCa(   0mS/cm^2),
               Kdr( 125mS/cm^2),
               H(   .05mS/cm^2),
               leak(.01mS/cm^2)];
PY_dynamics = HodgkinHuxley(Vₘ, PY_channels, gradients; geometry = geo);
@named PY = Compartment(LP_dynamics, extensions = [calcium_conversion]);

PY # display

topology = NetworkTopology([ABPD, LP, PY], [Glut, Chol]);
reversal_map = Dict([Glut => EGlut, Chol => EChol])

topology[ABPD, LP] = Glut(30nS)
topology[ABPD, LP] = Chol(30nS)
topology[ABPD, PY] = Glut(10nS)
topology[ABPD, PY] = Chol(3nS)
topology[LP, ABPD] = Glut(30nS)
topology[LP, PY] = Glut(1nS)
topology[PY, LP] = Glut(30nS)

network = NeuronalNetworkSystem(topology, reversal_map)
t_total = 10000
sim = Simulation(network, time = t_total*ms);
solution = solve(sim, RadauIIA5(), reltol=1e-5, abstol=1e-5);
# Plot at 5kHz sampling
fig = plot(solution(0.0:0.2:t_total; idxs = [ABPD.Vₘ, LP.Vₘ, PY.Vₘ]), size=(1200,800))
fig

# Uncomment and eval `png(...)` to save as PNG
# png(fig, "figure3e_simulated")

