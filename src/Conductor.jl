module Conductor

using DataStructures

using ModelingToolkit,
      Unitful,
      Unitful.DefaultSymbols,
      InteractiveUtils,
      IfElse,
      Symbolics,
      SymbolicUtils,
      Setfield,
      MacroTools

const MTK = ModelingToolkit

import Symbolics:
    get_variables,
    Symbolic,
    value,
    tosymbol,
    VariableDefaultValue,
    wrap, unwrap, Arr

import ModelingToolkit:
    toparam,
    isparameter,
    Equation,
    defaults,
    AbstractSystem,
    get_eqs,
    get_iv,
    get_ivs,
    get_states,
    get_observed,
    get_defaults,
    get_ps,
    get_systems,
    _merge,
    renamespace,
    hasdefault,
    getdefault,
    setdefault,
    AbstractTimeDependentSystem,
    independent_variables,
    get_variables!

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

export Gate, AlphaBeta, SteadyStateTau, SteadyState, ConstantValue,
       IonChannel, PassiveChannel, SynapticChannel, Synapse, Junction,
       AxialConductance

export EquilibriumPotential, Equilibrium, Equilibria, MembranePotential,
       IonCurrent, IonConcentration, Concentration, ExtrinsicPotential,
       Instrinsic, Extrinsic

export AuxConversion, D
export Simulation
export @named

export Calcium, Sodium, Potassium, Chloride, Cation, Anion, Leak, Ion, NonIonic

export t

export CompartmentSystem, Compartment,
       ConductanceSystem, Conductance, 
       NeuronalNetworkSystem, NeuronalNetwork,
       MultiCompartmentSystem, MultiCompartment

export output, get_output, timeconstant, steadystate, forward_rate,
       reverse_rate, hasexponent, exponent

export Sphere, Cylinder, Point, Unitless, area, radius, height


# Metadata IDs
struct ConductorUnits end # temporary shim until we implement MTK's unit checking
struct ConductorMaxConductance end

const ExprValues = Union{Expr,Symbol,Number}  # For use in macros
isfunction(ex::ExprValues) = try return eval(ex) isa Function catch; return false end

function extract_symbols(ex::ExprValues, out::Vector{Symbol}=[])
    if ~isfunction(ex) && isa(ex, Symbol)
        union!(out, [ex])
    end
    return ex
end

import Symbolics: unwrap, symtype, getindex_posthook

# Hijacked and modified from Symbolics.jl
function set_symarray_metadata(x, ctx, val)
    if symtype(x) <: AbstractArray
        if val isa AbstractArray
            getindex_posthook(x) do r,x,i...
                set_symarray_metadata(r, ctx, val[i...])
            end
        else
            getindex_posthook(x) do r,x,i...
                set_symarray_metadata(r, ctx, val)
            end
        end
    else
        setmetadata(x, ctx, val)
    end
end

include("primitives.jl")
include("gates.jl")
include("channels.jl")
include("compartments.jl")
include("multicompartment.jl")
include("networks.jl")
include("io.jl")
include("util.jl")

function Simulation(neuron::AbstractCompartmentSystem; time::Time, return_system = false)
    odesys = convert(ODESystem, neuron)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(odesys)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODAEProblem(simplified, [], (0., t_val), [])
    end
end

function Simulation(network::NeuronalNetworkSystem; time::Time, return_system = false)
    odesys = convert(ODESystem, network)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(odesys)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODAEProblem(simplified, [], (0., t_val), [])
    end
end

end # module
