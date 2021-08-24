# Gating variables (as an interface)

import ModelingToolkit:
    independent_variables

abstract type AbstractGatingSystem end # <: AbstractTimeDependentSystem end

struct GenericGate{S<:AbstractTimeDependentSystem} <: AbstractGatingSystem
    sys::S # symbolic system
    output::Num # the unitless output variable
    transform::Rule # output transformation on use (e.g. out^p)
end

struct Gate <: AbstractGatingSystem
    sys::ReactionSystem
    output::Num
    transform::Rule
end

Base.nameof(sys::AbstractGatingSystem) = nameof(getfield(sys, :sys))
ModelingToolkit.independent_variables(sys::AbstractGatingSystem) =
independent_variables(getfield(sys, :sys))
get_output(x::AbstractGatingSystem) = getfield(x, :output)
get_rule(x::AbstractGatingSystem) = getfield(x, :transform)
Base.merge(gate1::AbstractGatingSystem, gate2::AbstractGatingSystem) =
merge(getfield(gate1, :sys), getfield(gate2, :sys))

# TODO: hassteadystate(x::AbstractGatingSystem) = hasfield(typeof(x), :ss) ? !(isnothing(x.ss)) : false
# hasexponent(x::AbstractGatingSystem) = hasfield(typeof(x), :p) ? x.p !== one(Float64) : false

Base.convert(::Type{<:ODESystem}, x::Gate) = convert(ODESystem, x.sys; include_zero_odes=false)
# etc... for available conversions from Catalyst.jl

# Model types via trait types
abstract type AbstractGateModel end
struct SteadyStateTau <: AbstractGateModel end
struct AlphaBetaRates <: AbstractGateModel end

function Gate(::Type{AlphaBetaRates}, alpha::Num, beta::Num, p::Real;
        defaults::Dict=Dict(), symname::Symbol, name::Symbol = Base.gensym("Gate"))
    
    ss = alpha/(alpha + beta) # αₘ/(αₘ + βₘ)
    rule = @rule ~x => (~x)^p
    pars = Set{Num}()
    defaultmap = Dict()
    out = only(@variables $symname(t))
    states = Set{Num}([out])
    push!(defaultmap, out => ss) 

    rxns = [Reaction(alpha,nothing,[out]),
            Reaction(alpha+beta,[out],nothing)]

    for symbol in union(get_variables(alpha),get_variables(beta))
        isparameter(symbol) ? push!(pars, symbol) : push!(states, symbol)
        hasdefault(symbol) && push!(defaultmap, symbol => getdefault(symbol))
    end

    rxn_sys = ReactionSystem(rxns, t, states, pars; defaults = merge(defaultmap, defaults),
                             name=name)

    Gate(rxn_sys, out, rule)
end

function Gate(::Type{SteadyStateTau}, ss::Num, tau::Num, p::Real;
        defaults::Dict=Dict(), symname::Symbol, name::Symbol = Base.gensym("Gate"))
    
    rule = @rule ~x => (~x)^p
    pars = Set{Num}()
    defaultmap = Dict()
    out = only(@variables $symname(t))
    states = Set{Num}([out])
    push!(defaultmap, out => ss) 
    rxns = [Reaction(ss/tau,nothing,[out]),
            Reaction(1/tau,[out],nothing)]

    for symbol in union(get_variables(ss),get_variables(tau))
        isparameter(symbol) ? push!(pars, symbol) : push!(states, symbol)
        hasdefault(symbol) && push!(defaultmap, symbol => getdefault(symbol))
    end

    rxn_sys = ReactionSystem(rxns, t, states, pars; defaults = merge(defaultmap, defaults),
                             name=name)
    Gate(rxn_sys, out, rule)
end

function Gate(::Type{SteadyStateTau}; p = one(Int64), defaults::Dict=Dict(), kwargs...)
    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:m∞, :τₘ], syms)
        Gate(SteadyStateTau, kwargs[:m∞], kwargs[:τₘ], p; defaults=defaults, symname = :m)
    elseif issetequal([:h∞, :τₕ], syms)
        Gate(SteadyStateTau, kwargs[:h∞], kwargs[:τₕ], p; defaults=defaults, symname = :h)
    else
        throw("invalid keywords")
    end
end

function Gate(::Type{AlphaBetaRates}; p = one(Int64), defaults::Dict=Dict(), kwargs...)
    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:αₘ, :βₘ], syms)
        Gate(AlphaBetaRates, kwargs[:αₘ], kwargs[:βₘ], p; defaults=defaults, symname = :m)
    elseif issetequal([:αₕ, :βₕ], syms)
        Gate(AlphaBetaRates, kwargs[:αₕ], kwargs[:βₕ], p; defaults=defaults, symname = :h)
    elseif issetequal([:αₙ, :βₙ], syms)
        Gate(AlphaBetaRates, kwargs[:αₙ], kwargs[:βₙ], p; defaults=defaults, symname = :n)
    else
        throw("invalid keywords")
    end
end
