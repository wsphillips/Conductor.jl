# Gating variables (as an interface)

#abstract type AbstractKineticSystem end # <: AbstractTimeDependentSystem end
abstract type AbstractGatingVariable end

get_output(x::AbstractGatingVariable) = x.output
#has_transform(x::AbstractGatingSystem) = !(isnothing(getfield(x, :transform)))
#get_transform(x::AbstractGatingSystem) = getfield(x, :transform)
# TODO: hassteadystate(x::AbstractGatingSystem) = hasfield(typeof(x), :ss) ? !(isnothing(x.ss)) : false
# hasexponent(x::AbstractGatingSystem) = hasfield(typeof(x), :p) ? x.p !== one(Float64) : false

@enum GateVarType NullType Activation Inactivation
@enum GateVarForm AlphaBeta SteadyStateTau

abstract type AbstractGatingVariable

struct GatingVariable <: AbstractGatingVariable
    output::Num
    alpha::Num
    beta::Num
    steadystate::Num
    p::Real
    type::GateVarType # for future use
end

function GatingVariable(T::GateVarForm, x::Num, y::Num, p::Real = 1; name = Base.gensym("GateVar"))
    out = only(@variables $name(t))
    if T == AlphaBeta
        alpha, beta = x, y
        ss = alpha/(alpha + beta) 
    elseif T == SteadyStateTau
        ss, tau = x, y
        alpha = ss/tau         
        beta = inv(tau) - alpha
    end
    Gate(out, alpha, beta, ss, p, NullType)
end

# Likely replace this with macro
function GatingVariable(T::GateVarForm; p = one(Int64), kwargs...)
    syms = keys(kwargs)
    length(syms) !== 2 && throw("Invalid number of input equations.")
    if T == SteadyStateTau
        if issetequal([:m∞, :τₘ], syms)
            return Gate(SteadyStateTau, kwargs[:m∞], kwargs[:τₘ], p; name = :m)
        elseif issetequal([:h∞, :τₕ], syms)
            return Gate(SteadyStateTau, kwargs[:h∞], kwargs[:τₕ], p; name = :h)
        elseif isesetequal([:n∞, :τₙ], syms)
            return Gate(SteadyStateTau, kwargs[:n∞], kwargs[:τₙ], p; name = :n)
        else
            throw("invalid keywords")
        end
    elseif T == AlphaBeta
        if issetequal([:αₘ, :βₘ], syms)
            return Gate(AlphaBetaRates, kwargs[:αₘ], kwargs[:βₘ], p; name = :m)
        elseif issetequal([:αₕ, :βₕ], syms)
            return Gate(AlphaBetaRates, kwargs[:αₕ], kwargs[:βₕ], p; name = :h)
        elseif issetequal([:αₙ, :βₙ], syms)
            return Gate(AlphaBetaRates, kwargs[:αₙ], kwargs[:βₙ], p; name = :n)
        else
            throw("invalid keywords")
        end
    end
end

#=
struct GenericGate{S<:AbstractTimeDependentSystem} <: AbstractGatingSystem
    sys::S # symbolic system
    output::Num # the unitless output variable
    transform::Union{Nothing,Rule} # output transformation on use (e.g. out^p)
end
const Gate = GenericGate{ODESystem}
const ReactionGate = GenericGate{ReactionSystem}

struct Gate <: AbstractGatingSystem
    sys::ReactionSystem
    output::Num
    transform::Rule
end

#Base.nameof(sys::AbstractGatingSystem) = nameof(getfield(sys, :sys))
#independent_variables(sys::AbstractGatingSystem) = independent_variables(getfield(sys, :sys))
#Base.merge(gate1::AbstractGatingSystem, gate2::AbstractGatingSystem) =
#merge(getfield(gate1, :sys), getfield(gate2, :sys))
#Base.convert(::Type{<:ODESystem}, x::Gate) = convert(ODESystem, x.sys; include_zero_odes=false)
# etc... for available conversions from Catalyst.jl

=#

#=
function GatingVariable(::Type{AlphaBeta}, alpha::Num, beta::Num, p::Real; name::Symbol = Base.gensym("Gate"))
    
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

function GatingVariable(::Type{SteadyStateTau}, ss::Num, tau::Num, p::Real; name::Symbol = Base.gensym("Gate"))
    
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

function Gate(::Type{AlphaBeta}; p = one(Int64), defaults::Dict=Dict(), kwargs...)
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
=#
