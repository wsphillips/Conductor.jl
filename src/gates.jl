# Gating variables
abstract type AbstractGatingVariable end

output(x::AbstractGatingVariable) = getfield(x, :output)
timeconstant(x::AbstractGatingVariable) = getfield(x, :tau)
steadystate(x::AbstractGatingVariable) = getfield(x, :steadystate)
forward_rate(x::AbstractGatingVariable) = getfield(x, :alpha)
backward_rate(x::AbstractGatingVariable) = getfield(x: :beta)
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
    out = only(@variables $name(t))
    if T == AlphaBeta
        alpha, beta = x, y
        ss = alpha/(alpha + beta) 
        tau = inv(alpha + beta)
    elseif T == SteadyStateTau
        ss, tau = x, y
        alpha = ss/tau         
        beta = inv(tau) - alpha
    end
    GatingVariable(out, alpha, beta, ss, tau, p, NullType)
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


