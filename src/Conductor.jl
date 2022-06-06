    module Conductor

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

const ℱ = Unitful.q*Unitful.Na # Faraday's constant
const t = let name = :t; only(@variables $name) end
const D = Differential(t)
const ExprValues = Union{Expr,Symbol,Number}  # For use in macros

# Metadata IDs
struct ConductorUnits end # temporary shim until we implement MTK's unit checking
struct ConductorMaxConductance end

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

# Basic symbols
@enum PrimitiveSource Intrinsic Extrinsic

function MembranePotential(V0 = -60mV; dynamic = true, source::PrimitiveSource = Intrinsic,
                           n::Integer = 1, name::Symbol = :Vₘ)
    if isnothing(V0)
        if n == one(n) 
            ret = dynamic ? only(@variables $name(t)) :
                            only(@parameters $name)
        elseif n > one(n)
            ret = dynamic ? only(@variables $name[1:n](t)) :
                            only(@parameters $name[1:n])
        else
            throw("'n' must be greater than or equal to 1")
        end
    else
        V0_val = ustrip(Float64, mV, V0) #FIXME: assumes V0 <: Voltage
        if n == one(n)
            ret = dynamic ? only(@variables $name(t) = V0_val) :
                            only(@parameters $name = V0_val)
        elseif n > one(n)
            ret = dynamic ? only(@variables $name[1:n](t) = V0_val) :
                            only(@parameters $name[1:n] = V0_val)
        else
            throw("'n' must be greater than or equal to 1")
        end
    end

    ret = set_symarray_metadata(ret, PrimitiveSource, source)
    return ret
end

ExtrinsicPotential(;V0 = nothing, n = 1, dynamic = true,
                   source::PrimitiveSource = Extrinsic, name::Symbol = :Vₓ) = 
MembranePotential(V0; dynamic = dynamic, source = source, n = n, name = name)

@enum Location Outside Inside
# Custom Unitful.jl quantities
@derived_dimension SpecificConductance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^3 # conductance per unit area
@derived_dimension SpecificCapacitance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^4 # capacitance per unit area
@derived_dimension ConductancePerFarad 𝐓^-1 # S/F cancels out to 1/s; perhaps find a better abstract type?

#const TimeF64 = Quantity{Float64, 𝐓, U} where U

include("ions.jl")
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
