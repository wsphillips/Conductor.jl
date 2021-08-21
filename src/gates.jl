# Gating variables (as an interface)
abstract type AbstractGatingVariable end

# TODO: hassteadystate(x::AbstractGatingVariable) = hasfield(typeof(x), :ss) ? !(isnothing(x.ss)) : false
hasexponent(x::AbstractGatingVariable) = hasfield(typeof(x), :p) ? x.p !== one(Float64) : false
getsymbol(x::AbstractGatingVariable) = only(states(x.sys))
getequation(x::AbstractGatingVariable) = only(equations(x.sys))
getdefaults(x::AbstractGatingVariable) = defaults(x.sys)
subtype(x::AbstractGatingVariable) = typeof(x.sys)
subtype(x::Type{G}) where G<:AbstractGatingVariable = fieldtype(G,:sys)

struct Gate{S<:AbstractSystem} <: AbstractGatingVariable
    sys::S # symbolic system
    p::Float64 # optional exponent (defaults to 1)
end

abstract type AbstractGateModel end
struct SteadyStateTau <: AbstractGateModel end
struct AlphaBetaRates <: AbstractGateModel end

function Gate(::Type{SteadyStateTau}, name::Symbol, ss::Num, tau::Num, p::Float64;
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem) where S<:AbstractSystem
    sym = only(@variables $name(t))

    pars = Set{Num}()
    pars_defaults = Dict()
    for symbol in union(get_variables.([ss,tau])...)
        if isparameter(symbol)
            # Extra parameters are made unique to the gate ie. x => m₊x
            par = renamespace(name,symbol)
            # substitute in unique parameter
            alpha, beta = substitute.([alpha,beta], (symbol => par))
            # set parameter default
            haskey(defaults,symbol) ? push!(pars_defaults, par => defaults[symbol]) : nothing
            push!(pars,par)
        end
    end
    rn = ReactionSystem(
        [Reaction(ss/tau,nothing,[sym]),
        Reaction(1/tau,[sym],nothing)],
        t,[sym],pars;
        defaults=_merge(Dict(sym=>ss),pars_defaults),
        name=name
    )
    sys = convert(sys_type,rn)
    Gate(sys,p)
end

function Gate(::Type{AlphaBetaRates}, name::Symbol, alpha::Num, beta::Num, p::Float64;
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem) where S<:AbstractSystem
    sym = only(@variables $name(t))

    pars = Set{Num}()
    pars_defaults = Dict()
    for symbol in union(get_variables.([alpha,beta])...)
        if isparameter(symbol)
            # Extra parameters are made unique to the gate ie. x => m₊x
            par = renamespace(name,symbol)
            # substitute in unique parameter
            alpha, beta = substitute.([alpha,beta], (symbol => par))
            # set parameter default
            haskey(defaults,symbol) ? push!(pars_defaults, par => defaults[symbol]) : nothing
            push!(pars,par)
        end
    end
    ss = alpha/(alpha + beta) # αₘ/(αₘ + βₘ)
    rn = ReactionSystem(
        [Reaction(alpha,nothing,[sym]),
        Reaction(alpha+beta,[sym],nothing)],
        t,[sym],pars;
        defaults=_merge(Dict(sym=>ss),pars_defaults),
        name=name
    )
    sys = convert(sys_type,rn)
    Gate(sys,p)
end

Gate(t::Type{SteadyStateTau}, name::Symbol, alpha::Num, beta::Num, p::Real; kwargs...) =
Gate(t, name, alpha, beta, Float64(p); kwargs...)

Gate(t::Type{AlphaBetaRates}, name::Symbol, alpha::Num, beta::Num, p::Real; kwargs...) =
Gate(t, name, alpha, beta, Float64(p); kwargs...)

# TODO: find a nicer way to do this
function Gate(::Type{SteadyStateTau}; p = one(Float64),
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem, kwargs...) where S<:AbstractSystem

    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:m∞, :τₘ], syms)
        Gate(SteadyStateTau, :m, kwargs[:m∞], kwargs[:τₘ], p; defaults=defaults, sys_type=sys_type)
    elseif issetequal([:h∞, :τₕ], syms)
        Gate(SteadyStateTau, :h, kwargs[:h∞], kwargs[:τₕ], p; defaults=defaults, sys_type=sys_type)
    else
        throw("invalid keywords")
    end
end

function Gate(::Type{AlphaBetaRates}; p = one(Float64),
    defaults::Dict=Dict(), sys_type::Type{S}=ODESystem, kwargs...) where S<:AbstractSystem

    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:αₘ, :βₘ], syms)
        Gate(AlphaBetaRates, :m, kwargs[:αₘ], kwargs[:βₘ], p; defaults=defaults, sys_type=sys_type)
    elseif issetequal([:αₕ, :βₕ], syms)
        Gate(AlphaBetaRates, :h, kwargs[:αₕ], kwargs[:βₕ], p; defaults=defaults, sys_type=sys_type)
    elseif issetequal([:αₙ, :βₙ], syms)
        Gate(AlphaBetaRates, :n, kwargs[:αₙ], kwargs[:βₙ], p; defaults=defaults, sys_type=sys_type)
    else
        throw("invalid keywords")
    end
end
