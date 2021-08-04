# Classic Hodgkin Huxley neuron with a "current pulse" stimulus
using Conductor, IfElse, OrdinaryDiffEq, Unitful, ModelingToolkit
using Symbolics, ExprTools, MacroTools
using MacroTools: prettify, prewalk, rmlines, postwalk
using ExprTools: splitdef
using JuliaFormatter
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
@named simulation = ODESystem([D(neuron.sys.Isyn) ~ 0]; systems = [neuron.sys])
simp = structural_simplify(simulation)

(oop, inp) = Symbolics.build_function([x.rhs for x in equations(simp)],
                                   states(simp),
                                   ModelingToolkit.parameters(simp),
                                   independent_variable(simp))

cleanup = prettify(postwalk(simplify,prewalk(rmlines, inp)))

ivsuffix(x,sys) = endswith(string(x),"($(independent_variable(simp)))")

function chopiv(x,sys)
    suflen = length(string(independent_variable(sys))) + 2
    return Symbol(chop(string(x), tail = suflen))
end

cleanup = postwalk(x -> (x isa Symbol && ivsuffix(x,simp)) ? chopiv(x,simp) : x, cleanup)

def = ExprTools.splitdef(cleanup)
def[:name] = :hodgkin_huxley!

for (old, new) in zip(def[:args][1:3], [:du, :u, :p])
    cleanup = postwalk(x -> x isa Symbol && x == old ? new : x, cleanup)

#    def[:body] = postwalk(x -> x isa Symbol && x == old ? new : x, def[:body])
end

def[:args] = Any[:du,:u,:p,:t]

final = ExprTools.combinedef(def)

clipboard(format_text(string(final)))

hhfun = eval(final)
diff_vars = [true, false, true, true, true, true]


