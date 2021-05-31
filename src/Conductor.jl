module Conductor

using ModelingToolkit, Unitful, Unitful.DefaultSymbols, InteractiveUtils
using IfElse, Symbolics, SymbolicUtils

import Symbolics: get_variables, Symbolic, value, tosymbol, VariableDefaultValue
import ModelingToolkit: toparam, isparameter, Equation
import SymbolicUtils: FnType

import Unitful: Time, Voltage, Current, Molarity
import Unitful: mV, mS, cm, µF, mF, µm, pA, nA, mA, µA, ms, mM

import Base: show, display

export Gate, AlphaBetaRates, SteadyStateTau, IonChannel, PassiveChannel
export EquilibriumPotential, Equilibrium, Equilibria, MembranePotential, MembraneCurrent
export Soma, Simulation
export @named

const t = toparam(Num(Sym{Real}(:t)))
const D = Differential(t)

function symoft(name::Symbol)
    Num(Variable{FnType{Tuple{Any}, Real}}(name))(t)
end

const MembranePotential() = symoft(:Vₘ)

@derived_dimension SpecificConductance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^3 # conductance per unit area
@derived_dimension SpecificCapacitance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^4 # capacitance per unit area
@derived_dimension ConductancePerFarad 𝐓^-1 # S/F cancels out to 1/s; perhaps find a better abstract type?

# Ion species

@enum Location Outside Inside

abstract type Ion end
abstract type Cation <: Ion end
abstract type Anion <: Ion end

struct Calcium  <: Ion end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Chloride <: Ion end

# Convenience aliases
const Ca = Calcium
const Na = Sodium
const K = Potassium
const Cl = Chloride
const Mixed = Ion # non-specific ion
const Leak = Mixed

const PERIODIC_SYMBOL = IdDict(Na => :Na,
                               K  => :K,
                               Cl => :Cl,
                               Ca => :Ca,
                               Leak => :l)

export Calcium, Sodium, Potassium, Chloride, Cation, Anion, Leak, Ion

abstract type ConductorCurrentCtx end
abstract type ConductorEquilibriumCtx end
abstract type ConductorConcentrationCtx end

struct MembraneCurrent{I<:Ion,V<:Union{Nothing,Num,Symbolic,Current}}
    ion::Type{I}
    val::V
end

function MembraneCurrent{I}(val = nothing; name::Symbol = PERIODIC_SYMBOL[I]) where {I <: Ion}
    var = Sym{Real}(Symbol("I", name)) # FIXME: set this to symoft?
    var = setmetadata(var, ConductorCurrentCtx, MembraneCurrent(I, val))
    return val isa Current ? toparam(Num(var)) : Num(var)
end

abstract type AbstractConcentration end
abstract type AbstractIonGradient end

struct EquilibriumPotential{I<:Ion,V<:Union{Num,Symbolic,Voltage}} <: AbstractIonGradient
    ion::Type{I}
    val::V
end

const Equilibrium{I} = EquilibriumPotential{I} 

function EquilibriumPotential{I}(val; name::Symbol = PERIODIC_SYMBOL[I]) where {I <: Ion}
    var = Sym{Real}(Symbol("E", name))
    var = setmetadata(var, ConductorEquilibriumCtx, EquilibriumPotential(I, val))
    return val isa Voltage ? toparam(Num(var)) : Num(var)
end

struct IonConcentration{I<:Ion, L<:Location, V<:Union{Num,Symbolic,Molarity}} <:AbstractConcentration
    ion::Type{I}
    val::V
    loc::L
end

function Concentration(::Type{I}, val = 0mM, loc::Location = Inside, name::Symbol = PERIODIC_SYMBOL[I]) where {I <: Ion}
    var = Sym{Real}(Symbol("⟦",name,"⟧",(loc == Inside ? "ᵢ" : "ₒ")))
    var = setmetadata(var,  ConductorConcentrationCtx, IonConcentration(I, val, loc))
    return val isa Molarity ? toparam(Num(var)) : Num(var)
end

# Alternative constructor
function Equilibria(equil::Vector)
    out = Num[]
    for x in equil
        !(x.first <: Ion) && throw("Equilibrium potential must be associated with an ion type.")
        if typeof(x.second) <: Tuple
            tup = x.second
            typeof(tup[2]) !== Symbol && throw("Second tuple argument for $(x.first) must be a symbol.")
            push!(out, Equilibrium{x.first}(tup...))
        else
            push!(out, Equilibrium{x.first}(x.second))
        end
    end
    return out
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
    df::Equation # differential equation
    ss::Union{Nothing, Num} # optional steady-state expression for initialization
    p::Float64 # optional exponent (defaults to 1)
end

abstract type AbstractGateModel end
struct SteadyStateTau <: AbstractGateModel end
struct AlphaBetaRates <: AbstractGateModel end

function Gate(::Type{SteadyStateTau}, name::Symbol, ss::Num, tau::Num, p::Real)
    sym = symoft(name)
    df = D(sym) ~ (ss-sym)/tau # (m∞ - m)/τₘ
    return Gate(sym, df, ss, p)
end

function Gate(::Type{AlphaBetaRates}, name::Symbol, alpha::Num, beta::Num, p::Real)
    sym = symoft(name)
    df = D(sym) ~ alpha * (1 - sym) - beta*sym # αₘ(1 - m) - βₘ*m
    ss = alpha/(alpha + beta) # αₘ/(αₘ + βₘ)
    return Gate(sym, df, ss, p)
end

# TODO: find a nicer way to do this
function Gate(::Type{SteadyStateTau}; p = one(Float64), kwargs...)

    syms = keys(kwargs)

    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:m∞, :τₘ], syms)
        Gate(SteadyStateTau, :m, kwargs[:m∞], kwargs[:τₘ], p)
    elseif issetequal([:h∞, :τₕ], syms)
        Gate(SteadyStateTau, :h, kwargs[:h∞], kwargs[:τₕ], p)
    else
        throw("invalid keywords")
    end
end

function Gate(::Type{AlphaBetaRates}; p = one(Float64), kwargs...)

    syms = keys(kwargs)

    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:αₘ, :βₘ], syms)
        Gate(AlphaBetaRates, :m, kwargs[:αₘ], kwargs[:βₘ], p)
    elseif issetequal([:αₕ, :βₕ], syms)
        Gate(AlphaBetaRates, :h, kwargs[:αₕ], kwargs[:βₕ], p)
    elseif issetequal([:αₙ, :βₙ], syms)
        Gate(AlphaBetaRates, :n, kwargs[:αₙ], kwargs[:βₙ], p)
    else
        throw("invalid keywords")
    end
end

# Conductance types
abstract type AbstractConductance end

isbuilt(x::AbstractConductance) = x.sys !== nothing

mutable struct IonChannel <: AbstractConductance
    gbar::SpecificConductance # scaling term - maximal conductance per unit area
    conducts::DataType # ion permeability
    inputs::Vector{Num} # cell states dependencies (input args to kinetics); we can infer this
    params::Vector{Num}
    kinetics::Vector{<:AbstractGatingVariable} # gating functions; none = passive channel
    sys::Union{ODESystem, Nothing} # symbolic system
end

# Return ODESystem pretty printing for our wrapper types
Base.show(io::IO, ::MIME"text/plain", x::IonChannel) = Base.display(isbuilt(x) ? x.sys : x)

# General purpose constructor
function IonChannel(conducts::Type{I},
                    gate_vars::Vector{<:AbstractGatingVariable},
                    max_g::SpecificConductance;
                    name::Symbol, build = true) where {I <: Ion}
    
    # if no kinetics, the channel is just a scalar
    passive = length(gate_vars) == 0

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
    filter!(x -> !any(isequal(y, x) for y in gates), inputs)
    
    # g = total conductance (e.g. g(m,h) ~ ̄gm³h)
    if passive
        @variables g() # uncalled
    else
        @variables g(t)
    end

    states = [g
              gates
              inputs]

    params = @parameters gbar
    defaultmap = Pair[gbar => gbar_val]

    # the "output" of a channel is it's conductance: g
    if passive
        eqs = [g ~ gbar]
        push!(defaultmap, g => gbar_val)
    else
        geq = [g ~ gbar * prod(hasexponent(x) ? getsymbol(x)^x.p : getsymbol(x) for x in gate_vars)]
        eqs = [[getequation(x) for x in gate_vars]
               geq]
        for x in gate_vars
            push!(defaultmap, getsymbol(x) => hassteadystate(x) ? x.ss : 0.0) # fallback to zero
        end
    end
    system = build ? ODESystem(eqs, t, states, params; defaults = defaultmap, name = name) : nothing
    
    return IonChannel(max_g, conducts, inputs, params, gate_vars, system)
end

(chan::IonChannel)(newgbar::SpecificConductance) = (chan.gbar = newgbar)

# Alias for ion channel with static conductance
function PassiveChannel(conducts::Type{I}, max_g::SpecificConductance;
                        name::Symbol = Base.gensym(:Leak), build = true) where {I <: Ion}
    gate_vars = AbstractGatingVariable[] # ie. 'nothing'
    return IonChannel(conducts, gate_vars, max_g; name = name, build = build)
end

abstract type Geometry end
abstract type Sphere <: Geometry end
abstract type Cylinder <: Geometry end

struct Compartment{G}
    cap::SpecificCapacitance
    chans::Vector{<:AbstractConductance}
    states::Vector
    params::Vector
    sys::ODESystem
end

# TODO: input validation - check channels/reversals for duplicates/conflicts
function Compartment{Sphere}(channels::Vector{<:AbstractConductance},
                             gradients; name::Symbol, radius = 20µm,
                             capacitance::SpecificCapacitance = 1µF/cm^2,
                             V0::Voltage = -65mV,
                             applied::Current = 225nA,
                             #= auxstate = Equation[], =#
                             build = true)
   
    Vₘ = MembranePotential()
    r_val = ustrip(Float64, cm, radius)
    area = π*r_val^2
    
    params = @parameters cₘ aₘ Iapp()

    states = [Vₘ] # grow this as we discover/generate new states
    currents = [] 

    defaultmap = Pair[Iapp => ustrip(Float64, mA, applied),
                  aₘ => area,
                  Vₘ => ustrip(Float64, mV, V0),
                  cₘ => ustrip(Float64, mF/cm^2, capacitance)]

    eqs = Equation[]
    required_states = [] # states not produced or intrinsic (e.g. not currents or Vm)
                         # but that we discover referenced in the RHS
                         
    grad_meta = getmetadata.(gradients, ConductorEquilibriumCtx) 
    systems = []

    #= use for auxillary state transformations (e.g. net calcium current -> Ca concentration)
    for i in auxstate
        push!(states, i.lhs)
        append!(required_states, get_variables(i.rhs)...)
    end
    =# 

    # parse and build equations
    for chan in channels

        iontype = chan.conducts
        
        !isbuilt(chan) && build!(chan)

        sys = chan.sys
        push!(systems, sys) 
        
        # auto forward cell states to channels
        for inp in chan.inputs
            push!(required_states, inp)
            subinp = getproperty(sys, tosymbol(inp, escape=false))
            push!(eqs, inp ~ subinp)
            # Workaround for: https://github.com/SciML/ModelingToolkit.jl/issues/1013
            push!(defaultmap, subinp => inp)
        end
        
        # write the current equation state
        I = symoft(Symbol(:I,nameof(sys))) # alternatively "uncalled" term
        push!(states, I) 
        push!(currents, I)

        # for now, take the first reversal potential with matching ion type
        idx = findfirst(x -> x.ion == iontype, grad_meta)
        Erev = gradients[idx]
        eq = [I ~ aₘ * sys.g * (Vₘ - Erev)]
        rhs = grad_meta[idx].val 

        if typeof(rhs) <: Voltage
            push!(defaultmap, Erev => ustrip(Float64, mV, rhs))
            push!(params, Erev)
            push!(eqs, eq...) 
        else # if symbolic/dynamic reversal potentials (warning: not tested yet)
            push!(eq, Erev ~ rhs)
            rhs_vars = get_variables(rhs)
            filter!(x -> !isequal(x, value(Erev)), rhs_vars)
            append!(eqs, eq...)
            append!(required_states, rhs_vars)
        end
    end
    
    required_states = value.(required_states)
    unique!(required_states)

    # propagate default parameter values to channel systems
    vm_eq = D(Vₘ) ~ (Iapp - sum(currents))/(aₘ*cₘ)
    push!(eqs, vm_eq)
    system = build ? ODESystem(eqs, t, states, params; systems = systems, defaults = defaultmap, name = name) : nothing
    
    # Switch to outputting a compartment/neuron when we build-out networks/synapses
    return system#Compartment{Sphere}(capacitance, channels, states, params, system)
end

const Soma = Compartment{Sphere}

function Compartment{Cylinder}() end
# should also be able to parse "collections" of compartments" that have an adjacency list/matrix 

function Simulation(system; time::Time)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(system)
    return ODAEProblem(simplified, [], (0., t_val), [])
end

end # module

