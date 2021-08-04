# Classic Hodgkin Huxley neuron with a "current pulse" stimulus
using Conductor, IfElse, OrdinaryDiffEq, Unitful, ModelingToolkit, RuntimeGeneratedFunctions, JuliaFormatter
using Symbolics, ExprTools, MacroTools
using MacroTools: prettify, prewalk, rmlines, postwalk
using ExprTools: splitdef
#using JuliaFormatter
import Unitful: mV, mS, cm, µm, pA, nA, mA, µA, ms
import Conductor: Na, K # shorter aliases for Sodium/Potassium

Vₘ = MembranePotential()

nav_kinetics = [
    Gate(AlphaBetaRates,
         αₘ = IfElse.ifelse(Vₘ == -40.0, 1.0, (0.1*(Vₘ + 40.0))/(1.0 - exp(-(Vₘ + 40.0)/10.0))),
         βₘ = 4.0*exp(-(Vₘ + 65.0)/18.0),
         p = 3)
    Gate(AlphaBetaRates,
         αₕ = 0.07*exp(-(Vₘ+65.0)/20.0),
         βₕ = 1.0/(1.0 + exp(-(Vₘ + 35.0)/10.0)))]

kdr_kinetics = [
    Gate(AlphaBetaRates,
         αₙ = IfElse.ifelse(Vₘ == -55.0, 0.1, (0.01*(Vₘ + 55.0))/(1.0 - exp(-(Vₘ + 55.0)/10.0))),
         βₙ = 0.125 * exp(-(Vₘ + 65.0)/80.0),
         p = 4)]

@named NaV = IonChannel(Sodium, nav_kinetics, 120mS/cm^2) 
@named Kdr = IonChannel(Potassium, kdr_kinetics, 36mS/cm^2)
@named leak = PassiveChannel(Leak, 0.3mS/cm^2)

# Equilibrium potentials are a implicit description of a ion concentration gradient
gradients = Equilibria([Na   =>  50.0mV,
                        K    => -77.0mV,
                        Leak => -54.4mV])

area = 4*pi*(20µm)^2
pulse(t, current) = 100. < t < 200. ? ustrip(Float64, µA, 400pA) : 0.0
@register pulse(a,b)

@named neuron = Soma([NaV,Kdr,leak], gradients, stimulus = pulse, area = ustrip(Float64, cm^2, area));

t = 300 
sim = Simulation(neuron, time = t*ms)
############################################################################################

ivsuffix(x,sys) = endswith(string(x),"($(independent_variable(simp)))")

function chopiv(x,sys)
    suflen = length(string(independent_variable(sys))) + 2
    return Symbol(chop(string(x), tail = suflen))
end

function rgf_lambda_expr(x)
    args = first(ExprTools.parameters(typeof(x)))
    return Expr(:function, Expr(:tuple, args...), x.body)
end

rgf = sim.f.f
inpbody = prewalk(rmlines, rgf.body)
args = first(ExprTools.parameters(typeof(rgf)))

for (old, new) in zip(args, [:du, :u, :p, :t])
    inpbody = postwalk(x -> x isa Symbol && x == old ? new : x, inpbody)
end

inp = Expr(:function, :(hodgkin_huxley!(du, u, p, t)), inpbody)
inp = postwalk(x -> x isa RuntimeGeneratedFunction ? rgf_lambda_expr(x) : x, inp )
inp = prewalk(rmlines, inp)

cleanup = postwalk(x -> (x isa Symbol && ivsuffix(x,neuron.sys))?chopiv(x,neuron.sys):x,inp)

clipboard(format_text(string(prettify(cleanup))))


