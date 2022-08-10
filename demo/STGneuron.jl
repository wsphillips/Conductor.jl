# Multiple single neuron examples derived from the Prinz et al STG model.
# To try some of the other neuron variants, redefine the `channels` vector, using any of the
# additional commented out parameters listed below. Then re-run the neuron constructor,
# simulation, solve and plot as before.

using Conductor, OrdinaryDiffEq, Plots

include(joinpath(@__DIR__, "prinz_kinetics.jl"));

############################################################################################
# Prinz 2003 J. Neurophysiol
############################################################################################

# Figure 2C Bursting neuron
channels = [NaV(100mS/cm^2),
            CaT(0mS/cm^2),
            CaS(4mS/cm^2),
            KA(0mS/cm^2),
            KCa(15mS/cm^2),
            Kdr(50mS/cm^2),
            H(.02mS/cm^2),
            leak(.03mS/cm^2)]

dynamics = HodgkinHuxley(Vₘ, channels, gradients; geometry = geo);

@named neuron = CompartmentSystem(dynamics, extensions = [calcium_conversion]);
t_total = 5000.0
sim = Simulation(neuron, time = t_total*ms)
solution = solve(sim, Rosenbrock23(), abstol=1e-6, reltol=1e-6, saveat=0.2);

# Plot at 5kHz sampling
fig = plot(solution(0.0:0.5:t_total; idxs=[Vₘ]); size=(1200,800));
gui(fig)

# Uncomment and eval `png(...)` to save as PNG
# png(fig, "figure_simulated")

# Evaluate the alternate channel vectors below to experiment with different styles of
# neurons, derived directly from the original publications...

#= 
# Figure 3C Irregular bursting
channels = [NaV(400mS/cm^2),
            CaT(0mS/cm^2),
            CaS(8mS/cm^2),
            KA(50mS/cm^2),
            KCa(20mS/cm^2),
            Kdr(50mS/cm^2),
            H(.04mS/cm^2),
            leak(0mS/cm^2)]
=#

#=
# Figure 5A broad APs - qualitatively similar
channels = [NaV(200mS/cm^2),
            CaT(0mS/cm^2),
            CaS(2mS/cm^2),
            KA(0mS/cm^2),
            KCa(15mS/cm^2),
            Kdr(0mS/cm^2),
            H(.03mS/cm^2),
            leak(.04mS/cm^2)]
=#

#=
# Figure 5B "spike triplets" - doesn't reproduce; spikes repetitively
channels = [NaV(100mS/cm^2),
            CaT(0mS/cm^2),
            CaS(10mS/cm^2),
            KA(50mS/cm^2),
            KCa(10mS/cm^2),
            Kdr(50mS/cm^2),
            H(.03mS/cm^2),
            leak(.05mS/cm^2)]
=#

#=
# Figure 5C Plateau burster - doesnt reproduce, 
channels = [NaV(400mS/cm^2),
            CaT(2.5mS/cm^2),
            CaS(10mS/cm^2),
            KA(20mS/cm^2),
            KCa(5mS/cm^2),
            Kdr(25mS/cm^2),
            H(.04mS/cm^2),
            leak(.03mS/cm^2)]
=#

#=
# Figure 5D (Parabolic burster) - bursts but not sure if "parabolic"
channels = [NaV(100mS/cm^2),
            CaT(0mS/cm^2),
            CaS(6mS/cm^2),
            KA(10mS/cm^2),
            KCa(10mS/cm^2),
            Kdr(50mS/cm^2),
            H(.03mS/cm^2),
            leak(.05mS/cm^2)]
=#

#=
# figure 5G - Bursting with 17-15-17-15... APs - not sure if pattern matches
channels = [NaV(500mS/cm^2),
            CaT(2.5mS/cm^2),
            CaS(8mS/cm^2),
            KA(0mS/cm^2),
            KCa(15mS/cm^2),
            Kdr(75mS/cm^2),
            H(.05mS/cm^2),
            leak(0mS/cm^2)]
=#

############################################################################################
# Prinz 2003 - J. Neuroscience
############################################################################################

#=
# bursting neuron
channels = [NaV(200mS/cm^2),
            CaT(2.5mS/cm^2),
            CaS(4mS/cm^2),
            KA(50mS/cm^2),
            KCa(5mS/cm^2),
            Kdr(100mS/cm^2),
            H(.01mS/cm^2),
            leak(.01mS/cm^2)]
=#

#=
# tonic spiking neuron
channels = [NaV(200mS/cm^2),
            CaT(0mS/cm^2),
            CaS(4mS/cm^2),
            KA(10mS/cm^2),
            KCa(10mS/cm^2),
            Kdr(125mS/cm^2),
            H(.05mS/cm^2),
            leak(.04mS/cm^2)]
=#

############################################################################################
# Prinz 2004 - Nature Neurosci. -- (all below are intrinsically bursting ABPD)
############################################################################################

#=
# AB/PD 1
channels = [NaV(400mS/cm^2),
            CaT(2.5mS/cm^2),
            CaS(6mS/cm^2),
            KA(50mS/cm^2),
            KCa(10mS/cm^2),
            Kdr(100mS/cm^2),
            H(.01mS/cm^2),
            leak(0mS/cm^2)]
=#

#=
# AB/PD 2
channels = [NaV(100mS/cm^2),
            CaT(2.5mS/cm^2),
            CaS(6mS/cm^2),
            KA(50mS/cm^2),
            KCa(5mS/cm^2),
            Kdr(100mS/cm^2),
            H(.01mS/cm^2),
            leak(0mS/cm^2)]
=#

#=
# AB/PD 3
channels = [NaV(200mS/cm^2),
            CaT(2.5mS/cm^2),
            CaS(4mS/cm^2),
            KA(50mS/cm^2),
            KCa(5mS/cm^2),
            Kdr(50mS/cm^2),
            H(.01mS/cm^2),
            leak(0mS/cm^2)]
=#

#=
# AB/PD 4
channels = [NaV(200mS/cm^2),
            CaT(5mS/cm^2),
            CaS(4mS/cm^2),
            KA(40mS/cm^2),
            KCa(5mS/cm^2),
            Kdr(125mS/cm^2),
            H(.01mS/cm^2),
            leak(0mS/cm^2)]
=#

#=
# AB/PD 5
channels = [NaV(300mS/cm^2),
            CaT(2.5mS/cm^2),
            CaS(2mS/cm^2),
            KA(10mS/cm^2),
            KCa(5mS/cm^2),
            Kdr(125mS/cm^2),
            H(.01mS/cm^2),
            leak(0mS/cm^2)]
=#


