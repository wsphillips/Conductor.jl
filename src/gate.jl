# Gating variables
abstract type AbstractGatingVariable end

output(x::AbstractGatingVariable) = getfield(x, :output)

# Base overloads for indexing key-value pairs with dot syntax
function Base.getproperty(value::AbstractGatingVariable, name::Symbol)
    return getindex(getfield(value, :props), name)
end

function Base.setproperty!(value::AbstractGatingVariable, name::Symbol, x)
    return setindex!(getfield(value, :props), name, x)
end

function Base.propertynames(value::AbstractGatingVariable, private::Bool = false)
    props = collect(keys(getfield(value, :props))) 
    if private
        union!(props, fieldnames(typeof(value)))
    end
    return props
end

function Base.get(collection::AbstractGatingVariable, key::Symbol, default)
    return get(getfield(collection, :props), key, default)
end

"""
Abstract supertype for extending the behavior of `Gate`.

Stub subtypes of `GateVarForm` are used as traits when writing new methods that call `Gate`
objects.

# Example
```julia
struct MyNewGate <: GateVarForm end

function Conductor.get_eqs(g::Gate{MyNewGate}, comp::CompartmentSystem)
    ...
end
```
"""
abstract type GateVarForm end 

struct AlphaBeta <: GateVarForm end 
struct SteadyStateTau <: GateVarForm end
struct SteadyState <: GateVarForm end
struct ConstantValue <: GateVarForm end
struct HeavisideSum <: GateVarForm end

const SimpleGate = SteadyState

""" Test for Gate """
struct Gate{T<:GateVarForm} <: AbstractGatingVariable
    "something about output field"
    output::Num
    props::Dict{Symbol,Any}
end

function Gate{T}(output::Num; kwargs...) where T <: GateVarForm
    return Gate{T}(output, kwargs)
end

# Internal API: Gate property getters
timeconstant(x::Gate{<:Union{AlphaBeta,SteadyStateTau}}) = x.tau
steadystate(x::Gate) = x.ss
forward_rate(x::Gate{<:Union{AlphaBeta,SteadyStateTau}}) = x.alpha
reverse_rate(x::Gate{<:Union{AlphaBeta,SteadyStateTau}}) = x.beta
Base.exponent(x::AbstractGatingVariable) = get(x, :p, 1)

"""
    Gate(::Type{AlphaBeta}, alpha, beta; name = Base.gensym("GateVar"))

A gate that accepts expressions for forward (α) and reverse (β) reaction rates as
descriptors for its kinetics.

See also: [`get_eqs`](@ref).
"""
function Gate(::Type{AlphaBeta}, alpha, beta; name = Base.gensym("GateVar"), kwargs...)
    ss = alpha/(alpha + beta)
    tau = inv(alpha + beta)
    out = only(@variables $name(t) = ss)
    Gate{AlphaBeta}(out; alpha = alpha, beta = beta, ss = ss, tau = tau, kwargs...)
end

"""
    Gate(::Type{SteadyStateTau}, ss, tau; name = Base.gensym("GateVar"))

A gate that accepts expressions for its steady-state activation, x∞(Vₘ), and the time
constant, τₓ(Vₘ), as descriptors for its kinetics.

See also: [`get_eqs`](@ref).
"""
function Gate(::Type{SteadyStateTau}, ss, tau; name = Base.gensym("GateVar"), kwargs...)
    alpha = ss/tau
    beta = inv(tau) - alpha
    out = only(@variables $name(t) = ss)
    return Gate{SteadyStateTau}(out; alpha = alpha, beta = beta, ss = ss, tau = tau, kwargs...)
end

function Base.convert(::Type{Gate{SteadyState}},
                      x::Union{Gate{AlphaBeta},Gate{SteadyStateTau}})
    return Gate(SteadyState, steadystate(x), p = exponent(x),
                name = Symbolics.tosymbol(output(x), escape=false))
end

"""
    Gate(::Type{SimpleGate}, ss; name = Base.gensym("GateVar"))


"""
function Gate(::Type{SteadyState}, ss; name = Base.gensym("GateVar"), kwargs...)
    out = only(@variables $name(t) = ss)
    return Gate{SteadyState}(out; ss = ss, kwargs...)
end

"""
    Gate(::Type{ConstantValue}, val; name = Base.gensym("GateVar"))


"""
function Gate(::Type{ConstantValue}, val; name = Base.gensym("GateVar"), kwargs...)
    out = only(@parameters $name = val)
    return Gate{ConstantValue}(out; val = val, kwargs...)
end

"""

"""
function Gate(::Type{HeavisideSum}, threshold = 0mV, saturation = 125;
              name = Base.gensym("GateVar"), kwargs...) 
    out = only(@variables $name(t) = 0.0) # synaptically activated gate inits to 0.0
    return Gate{HeavisideSum}(out; threshold = threshold, saturation = saturation,
                              kwargs...)
end 

"""
    get_eqs(var::Gate{HeavisideSum}, chan)


"""
function get_eqs(var::Gate{HeavisideSum}, chan)
    thold, sat = var.threshold, var.saturation
    thold_val = ustrip(Float64, mV, thold)
    out = output(var)
    isempty(subscriptions(chan)) && return [D(out) ~ 0]
    @named Vₓ = ExtrinsicPotential(n = length(subscriptions(chan))) 
    # Derived from Pinsky & Rinzel 1994 - Equation 4 
    # S'ᵢ = ∑ 𝐻(Vⱼ - 10) - Sᵢ/150
    return[D(out) .~ sum(ModelingToolkit.scalarize(Vₓ .>= thold_val) .- (out/sat))]
end

"""
    get_eqs(var::Gate{<:Union{AlphaBeta, SteadySateTau}}, chan)

Generate the voltage- and time-dependent differential equation modeling the output of a gate
that was specified with [α(Vₘ),β(Vₘ)] or [x∞(Vₘ),τₓ(Vₘ)].

For 
"""
function get_eqs(var::Gate{<:Union{AlphaBeta,SteadyStateTau}}, chan)
    x, x∞, τₓ = output(var), steadystate(var), timeconstant(var)
    return [D(x) ~ inv(τₓ)*(x∞ - x)]
end

"""
    get_eqs(var::Gate{SimpleGate}, chan)


"""
get_eqs(var::Gate{SteadyState}, chan) = [output(var) ~ steadystate(var)]

"""
    get_eqs(var::Gate{ConstantValue}, chan)


"""
get_eqs(var::Gate{ConstantValue}, chan) = Equation[]

############################################################################################
# Macros (needs updating)
############################################################################################

macro gate(ex::Expr, p::Expr = :(p=1))
    name = Base.gensym("GateVar")
    _make_gate_variable(name, ex, p)
end

macro gate(name::Symbol, ex::Expr, p::Expr = :(p=1))
    _make_gate_variable(name, ex, p)
end

function _make_gate_variable(name::Symbol, ex::Expr, p::Expr = :(p=1))
    ex = MacroTools.striplines(ex)
    length(ex.args) !== 2 && throw("Invalid number of input equations.")
    p.args[1] != :p && throw("Please use `p` to define a gate exponent.")
    gate_pow = p.args[2]

    rate_vars = map(eq -> eq.args[1], ex.args)

    #= Find paramters in rate equations =#
    param_list = Symbol[]
    foreach(rate->MacroTools.postwalk(x -> extract_symbols(x, param_list), rate.args[2]), ex.args)
    filter!(x -> x != :Vₘ, param_list)

    if issetequal([:α, :β], rate_vars)  # Equation order doesn't matter
        rate_type = AlphaBeta
        gate_expr = :(Gate($AlphaBeta, α, β, $gate_pow; name=$(QuoteNode(name))))
    elseif issetequal([:ss, :τ], rate_vars)  # Equation order doesn't matter
        rate_type = SteadyStateTau
        gate_expr = :(Gate($SteadyStateTau, ss, τ, $p.args[2]; name=$(QuoteNode(name))))
    else
        throw("invalid keywords")
    end
    pushfirst!(ex.args, :(Vₘ = MembranePotential()))
    foreach(param -> pushfirst!(ex.args, :(@parameters $param)), param_list)  # generate parameters
    push!(ex.args, gate_expr)
    ex
end

