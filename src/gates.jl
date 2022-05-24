# Gating variables
abstract type AbstractGatingVariable end

output(x::AbstractGatingVariable) = getfield(x, :output)
timeconstant(x::AbstractGatingVariable) = getfield(x, :tau)
steadystate(x::AbstractGatingVariable) = getfield(x, :steadystate)
forward_rate(x::AbstractGatingVariable) = getfield(x, :alpha)
reverse_rate(x::AbstractGatingVariable) = getfield(x, :beta)
hasexponent(x::AbstractGatingVariable) = hasfield(typeof(x), :p) ? getfield(x, :p) !== one(typeof(x.p)) : false
exponent(x::AbstractGatingVariable) = getfield(x, :p)

abstract type GateVarForm end 
struct AlphaBeta <: GateVarForm end 
struct SteadyStateTau <: GateVarForm end
struct SteadyState <: GateVarForm end
struct ConstantValue <: GateVarForm end

struct Gate{T<:GateVarForm} <: AbstractGatingVariable
    output::Num
    alpha::Num
    beta::Num
    steadystate::Num
    tau::Num
    p::Real
end

function Gate(::Type{AlphaBeta}, x, y, p = 1; name = Base.gensym("GateVar"))
    alpha, beta = x, y
    ss = alpha/(alpha + beta)
    tau = inv(alpha + beta)
    out = only(@variables $name(t) = ss)
    Gate{AlphaBeta}(out, alpha, beta, ss, tau, p)
end

function Gate(::Type{SteadyStateTau}, x, y, p = 1; name = Base.gensym("GateVar"))
    ss, tau = x, y
    alpha = ss/tau
    beta = inv(tau) - alpha
    out = only(@variables $name(t) = ss)
    return Gate{SteadyStateTau}(out, alpha, beta, ss, tau, p)
end

function Base.convert(::Type{Gate{SteadyState}}, x::Union{Gate{AlphaBeta},Gate{SteadyStateTau}})
    return Gate(SteadyState, steadystate(x), exponent(x), name = Symbolics.tosymbol(output(x), escape=false))
end

function Gate(::Type{SteadyState}, x, p = 1; name = Base.gensym("GateVar"))
    out = only(@variables $name(t) = x)
    return Gate{SteadyState}(out, 0, 0, x, 0, p)
end

function Gate(::Type{ConstantValue}, x, p = 1; name = Base.gensym("GateVar"))
    out = only(@parameters $name = x)
    return Gate{ConstantValue}(out, 0, 0, x, 0, p)
end

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

