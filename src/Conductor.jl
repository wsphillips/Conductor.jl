module Conductor

using ModelingToolkit,
      Catalyst,
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
    wrap

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
    AbstractSystem,
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

export Gate, AlphaBeta, SteadyStateTau, IonChannel, PassiveChannel, SynapticChannel, Synapse
export EquilibriumPotential, Equilibrium, Equilibria, MembranePotential, IonCurrent
export AuxConversion, D, NeuronalNetwork
export Simulation, Concentration, IonConcentration
export @named
export Calcium, Sodium, Potassium, Chloride, Cation, Anion, Leak, Ion, NonIonic
export t
export CompartmentSystem, ConductanceSystem, NeuronalNetworkSystem, Conductance, Compartment
export output, get_output, timeconstant, steadystate, forward_rate, reverse_rate, hasexponent, exponent
export Sphere, Cylinder, Point, area, radius, height

const ℱ = Unitful.q*Unitful.Na # Faraday's constant
const t = let name = :t; only(@parameters $name) end
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

# Basic symbols
function MembranePotential(V0 = -60mV; dynamic = true, name::Symbol = :Vₘ)
    # remove when we have proper unit checks
    V0_val = ustrip(Float64, mV, V0)
    return dynamic ? only(@variables $name(t) = V0_val) : only(@parameters $name = V0_val)
end

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
include("networks.jl")
include("io.jl")

#=
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
=#

function Simulation(neuron::CompartmentSystem; time::Time, return_system = false)
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
