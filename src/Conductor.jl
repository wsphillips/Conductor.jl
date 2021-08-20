module Conductor

using ModelingToolkit,       
      Unitful,               
      Unitful.DefaultSymbols,
      InteractiveUtils,
      IfElse,                
      Symbolics,             
      SymbolicUtils,         
      Setfield

import Symbolics:
    get_variables,
    Symbolic,
    value,
    tosymbol,
    VariableDefaultValue,
    wrap

import ModelingToolkit:
    toparam,
    isparameter,
    Equation,
    defaults,
    AbstractSystem,
    get_eqs,
    get_states,
    get_observed,
    get_defaults,
    get_ps,
    get_systems

import ModelingToolkit.SciMLBase: parameterless_type

import Unitful:
    Time,
    Voltage,
    Current,
    Molarity,
    ElectricalConductance

import SymbolicUtils: FnType
import Unitful: mV, mS, cm, µF, mF, µm, pA, nA, mA, µA, ms, mM, µM
import Base: show, display

export Gate, AlphaBetaRates, SteadyStateTau, IonChannel, PassiveChannel, SynapticChannel
export EquilibriumPotential, Equilibrium, Equilibria, MembranePotential, MembraneCurrent
export AuxConversion, D, Network
export Soma, Simulation, Concentration, IonConcentration
export @named
export Calcium, Sodium, Potassium, Chloride, Cation, Anion, Leak, Ion

const ℱ = Unitful.q*Unitful.Na # Faraday's constant
const t = let name = :t; only(@parameters $name) end
const D = Differential(t)

# Temporary fix until https://github.com/SciML/ModelingToolkit.jl/issues/1223
# gets resolved. Non-flattening form of extend
function _extend(sys::AbstractSystem, basesys::AbstractSystem; name::Symbol=nameof(sys))
    T = parameterless_type(basesys)
    ivs = independent_variables(basesys)
    if !(typeof(sys) <: T)
        if length(ivs) == 0
            sys = convert_system(T, sys)
        elseif length(ivs) == 1
            sys = convert_system(T, sys, ivs[1])
        else
            throw("Extending multivariate systems is not supported")
        end
    end

    eqs = union(get_eqs(basesys), get_eqs(sys))
    sts = union(get_states(basesys), get_states(sys))
    ps = union(get_ps(basesys), get_ps(sys))
    obs = union(get_observed(basesys), get_observed(sys))
    defs = merge(get_defaults(basesys), get_defaults(sys)) # prefer `sys`
    syss = union(get_systems(basesys), get_systems(sys))

    if length(ivs) == 0
        T(eqs, sts, ps, observed = obs, defaults = defs, name=name, systems = syss)
    elseif length(ivs) == 1
        T(eqs, ivs[1], sts, ps, observed = obs, defaults = defs, name = name, systems = syss)
    end
end

# Helper utils -- FIXME: MTK has these now
#hasdefault(x::Symbolic) = hasmetadata(x, VariableDefaultValue) ? true : false
#hasdefault(x::Num) = hasdefault(ModelingToolkit.value(x))
#hasdefault(x) = false    
#
#getdefault(x::Symbolic) = hasdefault(x) ? getmetadata(x, VariableDefaultValue) : nothing
#getdefault(x::Num) = getdefault(ModelingToolkit.value(x))

# Basic symbols
function MembranePotential()
    name = :Vₘ
    return only(@variables $name(t))
end

@enum Location Outside Inside

# Custom Unitful.jl quantities
@derived_dimension SpecificConductance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^3 # conductance per unit area
@derived_dimension SpecificCapacitance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^4 # capacitance per unit area
@derived_dimension ConductancePerFarad 𝐓^-1 # S/F cancels out to 1/s; perhaps find a better abstract type?

include("ions.jl")
include("gates.jl")
include("channels.jl")
include("compartments.jl")
include("networks.jl")
include("io.jl")

function Simulation(network; time::Time, system = false)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(network)
    if system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODAEProblem(simplified, [], (0., t_val), [])
    end
end

function Simulation(neuron::Soma; time::Time, system = false)
    # for a single neuron, we just need a thin layer to set synaptic current constant
    old_sys = neuron.sys
    Isyn = getproperty(old_sys, :Isyn, namespace=false)
    wrapper = ODESystem([D(Isyn) ~ 0]; name = Base.gensym(:wrapper))
    # Use non-flattening extend
    new_neuron_sys = _extend(wrapper, old_sys; name = nameof(old_sys))
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(new_neuron_sys)
    if system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODAEProblem(simplified, [], (0., t_val), [])
    end
end

end # module
