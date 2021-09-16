# Gating variables
abstract type AbstractGatingVariable end

output(x::AbstractGatingVariable) = getfield(x, :output)
timeconstant(x::AbstractGatingVariable) = getfield(x, :tau)
steadystate(x::AbstractGatingVariable) = getfield(x, :steadystate)
forward_rate(x::AbstractGatingVariable) = getfield(x, :alpha)
reverse_rate(x::AbstractGatingVariable) = getfield(x, :beta)
hasexponent(x::AbstractGatingVariable) = hasfield(typeof(x), :p) ? getfield(x, :p) !== one(typeof(x.p)) : false
exponent(x::AbstractGatingVariable) = hasexponent(x) ? getfield(x, :p) : nothing

@enum GateVarType NullType Activation Inactivation
@enum GateVarForm AlphaBeta SteadyStateTau

struct GatingVariable <: AbstractGatingVariable
    output::Num
    alpha::Num
    beta::Num
    steadystate::Num
    tau::Num
    p::Real
    type::GateVarType # for future use
end

const Gate = GatingVariable

function GatingVariable(T::GateVarForm, x::Num, y::Num, p::Real = 1; name = Base.gensym("GateVar"))
    if T == AlphaBeta
        alpha, beta = x, y
        ss = alpha/(alpha + beta)
        tau = inv(alpha + beta)
    elseif T == SteadyStateTau
        ss, tau = x, y
        alpha = ss/tau
        beta = inv(tau) - alpha
    end

    out = only(@variables $name(t) = ss)
    GatingVariable(out, alpha, beta, ss, tau, p, NullType)
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
        gate_expr = :(GatingVariable($AlphaBeta, α, β, $gate_pow; name=$(QuoteNode(name))))
    elseif issetequal([:ss, :τ], rate_vars)  # Equation order doesn't matter
        rate_type = SteadyStateTau
        gate_expr = :(GatingVariable($SteadyStateTau, ss, τ, $p.args[2]; name=$(QuoteNode(name))))
    else
        throw("invalid keywords")
    end
    pushfirst!(ex.args, :(Vₘ = MembranePotential()))
    foreach(param -> pushfirst!(ex.args, :(@parameters $param)), param_list)  # generate parameters
    push!(ex.args, gate_expr)
    ex
end

# Likely replace this with macro
# function GatingVariable(T::GateVarForm; p = one(Int64), kwargs...)
#     syms = keys(kwargs)
#     length(syms) !== 2 && throw("Invalid number of input equations.")
#     if T == SteadyStateTau
#         if issetequal([:m∞, :τₘ], syms)
#             return Gate(SteadyStateTau, kwargs[:m∞], kwargs[:τₘ], p; name = :m)
#         elseif issetequal([:h∞, :τₕ], syms)
#             return Gate(SteadyStateTau, kwargs[:h∞], kwargs[:τₕ], p; name = :h)
#         elseif isesetequal([:n∞, :τₙ], syms)
#             return Gate(SteadyStateTau, kwargs[:n∞], kwargs[:τₙ], p; name = :n)
#         else
#             throw("invalid keywords")
#         end
#     elseif T == AlphaBeta
#         if issetequal([:αₘ, :βₘ], syms)
#             return Gate(AlphaBeta, kwargs[:αₘ], kwargs[:βₘ], p; name = :m)
#         elseif issetequal([:αₕ, :βₕ], syms)
#             return Gate(AlphaBeta, kwargs[:αₕ], kwargs[:βₕ], p; name = :h)
#         elseif issetequal([:αₙ, :βₙ], syms)
#             return Gate(AlphaBeta, kwargs[:αₙ], kwargs[:βₙ], p; name = :n)
#         else
#             throw("invalid keywords")
#         end
#     end
# end
