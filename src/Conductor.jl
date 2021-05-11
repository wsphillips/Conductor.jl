module Conductor

using ModelingToolkit, Unitful, Unitful.DefaultSymbols, InteractiveUtils, Symbolics
import Symbolics: get_variables, Symbolic, value
import Unitful: Voltage, mV, mS, cm, µF, mF
import Base: show, display

import ModelingToolkit: toparam
import ModelingToolkit.SymbolicUtils: FnType

const t = toparam(Num(Sym{Real}(:t)))
const D = Differential(t)

function symoft(name::Symbol)
    Num(Variable{FnType{Tuple{Any}, Real}}(name))(t)
end

const MembranePotential = symoft(:Vₘ)

@derived_dimension SpecificConductance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^3 # conductance per unit area
@derived_dimension SpecificCapacitance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^4 # capacitance per unit area
@derived_dimension ConductancePerFarad 𝐓^-1 # S/F cancels out to 1/s; perhaps find a better abstract type?

# Ion species
abstract type Ion end
struct Calcium <: Ion end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Chloride <: Ion end
struct Mixed <: Ion end # placeholder for non-specific ion
export Calcium, Sodium, Potassium, Chloride, Mixed

# Convenience aliases
const Ca = Calcium()
const Na = Sodium()
const K = Potassium()
const Cl = Chloride()

# These aren't used yet--not settled on how best to do this...
struct MembraneCurrent{I<:Ion} <: Real
    val::Union{Num,Symbolic}
end

struct EquilibriumPotential{I<:Ion} <: Real
    val::Union{Num, Symbolic}
end

"""
AbstractGatingVariable, in the most generic case, is any function that returns a
dimensionless value (weight).
"""
abstract type AbstractGatingVariable end

hassteadystate(x::AbstractGatingVariable) = hasfield(typeof(x), :ss) ? !(isnothing(x.ss)) : false
hasexponent(x::AbstractGatingVariable) = hasfield(typeof(x), :p) ? x.p !== one(Float64) : false
getsymbol(x::AbstractGatingVariable) = x.sym
getequation(x::AbstractGatingVariable) = x.df

struct Gate <: AbstractGatingVariable
    sym::Num # symbol/name (e.g. m, h)
    df::Symbolics.Equation # differential equation
    ss::Union{Nothing, Num} # optional steady-state expression for initialization
    p::Float64 # optional exponent (defaults to 1)
end

abstract type AbstractGateModel end
struct SteadyStateTau <: AbstractGateModel end
struct AlphaBetaRates <: AbstractGateModel end

function Gate(::SteadyStateTau, name::Symbol, ss::Num, tau::Num, p::Real)
    sym = symoft(name)
    df = D(sym) ~ (ss-sym)/tau # (m∞ - m)/τₘ
    return Gate(sym, df, ss, p)
end

function Gate(::AlphaBetaRates, name::Symbol, alpha::Num, beta::Num, p::Real)
    sym = symoft(name)
    df = D(sym) ~ alpha * (1 - sym) - beta*sym # αₘ(1 - m) - βₘ*m
    ss = alpha/(alpha + beta) # αₘ/(αₘ + βₘ)
    return Gate(sym, df, ss, p)
end

# TODO: find a nicer way to do this?
function Gate(::SteadyStateTau; p = one(Float64), kwargs...)

    syms = keys(kwargs)

    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:m∞, :τₘ], syms)
        Gate(SteadyStateTau(), :m, kwargs[:m∞], kwargs[:τₘ], p)
    elseif issetequal([:h∞, :τₕ], syms)
        Gate(SteadyStateTau(), :h, kwargs[:h∞], kwargs[:τₕ], p)
    else
        throw("unrecognized keywords")
    end
end

function Gate(::AlphaBetaRates; p = one(Float64), kwargs...)

    syms = keys(kwargs)

    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:αₘ, :βₘ], syms)
        Gate(SteadyStateTau(), :m, kwargs[:αₘ], kwargs[:βₘ], p)
    elseif issetequal([:αₕ, :βₕ], syms)
        Gate(SteadyStateTau(), :h, kwargs[:αₕ], kwargs[:βₕ], p)
    else
        throw("unrecognized keywords")
    end
end

# Conductance types
abstract type AbstractConductance end

mutable struct IonChannel <: AbstractConductance
    gbar::SpecificConductance # scaling term - maximal conductance per unit area
    conducts::Vector{<:Ion} # ion permeability
    inputs::Vector{Num} # cell states dependencies (input args to kinetics); we can infer this
    params::Vector{Num}
    kinetics::Vector{<:AbstractGatingVariable} # gating functions; none = passive channel
    sys::Union{ODESystem, Nothing} # symbolic system
end

# General purpose constructor
function IonChannel(conducts::Vector{<:Ion},
                    gate_vars::Vector{<:AbstractGatingVariable},
                    max_g::SpecificConductance;
                    name::Symbol, build = true)
    
    # TODO: Generalize to other possible units (e.g. S/F)
    gbar_val = ustrip(Float64, mS/cm^2, max_g)
    
    gates = [getsymbol(x) for x in gate_vars]

    # retrieve all variables present in the RHS of kinetics equations
    inputs = []

    for i in gate_vars
        syms = get_variables(getequation(i))
        for j in syms
            push!(inputs, j)
        end
    end
    
    # filter duplicates + self references from results
    unique!(inputs)
    inputs = filter(x -> !any(isequal(y, x) for y in gates), inputs)

    states = [@variables(g) # g = total conductance (e.g. ̄gm³h)
              gates
              inputs]

    params = @parameters gbar=gbar_val
    
    # the "output" of a channel is it's conductance: g
    eqs = [[getequation(x) for x in gate_vars]
           [g ~ gbar * prod(hasexponent(x) ? getsymbol(x)^x.p : getsymbol(x) for x in gate_vars)]]

    system = build ? ODESystem(eqs, t, states, params; name = name) : nothing
    
    return IonChannel(max_g, conducts, inputs, params, gate_vars, system)
end

# Return ODESystem pretty printing for our wrapper types
Base.show(io::IO, ::MIME"text/plain", x::IonChannel) = Base.display(isbuilt(x) ? x.sys : x)

IonChannel(ion::I, args...; kwargs...) where {I <: Ion} = IonChannel([ion], args...; kwargs...)

# Alias for ion channel with static conductance
function PassiveChannel(
    ions::Vector{<:Ion},
    max_g::SpecificConductance;
    name::Symbol = Base.gensym("Leak"),
    build = false)

    # Strip off units
    gbar_val = ustrip(Float64, mS/cm^2, max_g)

    states = @variables g()
    params = @parameters gbar
    
    eqs = [g ~ gbar]

    defaultmap = Dict(gbar => gbar_val)
    system = build ? ODESystem(eqs, t, states, params; name = name, defaults = defaultmap) : nothing

    return IonChannel(max_g, ions, [], [], defaultmap, system)
end

isbuilt(x::AbstractConductance) = x.sys !== nothing

# flags to know: VariableDefaultValue, MTKParameterCtx 
# Symbolics.setmetadata(x::Symbolic, t::DataType, value)
# Symbolics.makesym(MembranePotential, escape=false)
#=
struct InfTau <: AbstractGatingVariable
    sym::Num
    g∞::Function
    τg::Function
    p::Real # exponent
end

struct AlphaBeta <: AbstractGatingVariable
    sym::Num
    α::Function
    β::Function
    p::Real # exponent
end

AlphaBeta(α, β; exponent = 1) = AlphaBeta(α, β, exponent)
SSTau(g∞, τg; exponent = 1) = SSTau(g∞, τg, exponent)

# Steady-state channel gating
steady_state(state::AlphaBeta, V) = state.α(V)/(state.α(V) + state.β(V))
steady_state(state::InfTau, V) = state.g∞(V)

# Differential equations
function gate_equation(model::AlphaBeta, gate::Num, V::Num)
    α, β = model.α, model.β
    D(gate) ~ α(V) * (1 - gate) - β(V)*gate
end

function gate_equation(model::InfTau, gate::Num, V::Num)
    τ, g∞ = model.τg, model.g∞
    D(gate) ~ (1/τ(V))*(g∞(V) - gate)
end

# Optionally make gating variables callable
#(model::AlphaBeta)(gate::Num, V::Num) = gate_equation(model, gate, V)
#(model::SSTau)(gate::Num, V::Num) = gate_equation(model, gate, V)
=#

#=
abstract type AbstractNeuron end

struct Neuron <: AbstractNeuron
    sys::ODESystem
end

function Neuron(channels::Vector{<:AbstractConductance}, reversals; cm::SpecificCapacitance = 1µF/cm^2) 
    
    #Strip off units
    area = 1.0 # for now assume surface area 1.0
    cm_val = ustrip(Float64, mF/cm^2, cm)
    rev_vals = IdDict(reversals) 
    states = @variables V(t), I[1:length(channels)](t)
    params = @parameters Cm E[1:length(reversals)]
    
    ions = unique([ion_type(x) for x in channels])
    eq = [D(V) ~ -sum(chan.sys.g * area * (V - E[]))]

end
=#
# Helper because we can't use ∉(x::Num, itr) but isequal works...
function remove_symbol_set(badset, collection)
    newset = []
    for x in collection
        good = true
        for y in badset
            isequal(x,y) && (good = false) 
        end
        good == true && push!(newset, x)
    end
    return newset
end
=#

end # module
