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
$(TYPEDEF)

Abstract supertype for extending the behavior of `Gate`.

Stub subtypes of `GateVarForm` are used as traits when writing new methods that call `Gate`
objects.

# Example
```julia
struct MyNewGate <: Conductor.GateVarForm end

function Conductor.get_eqs(g::Gate{MyNewGate}, comp::CompartmentSystem)
    # a function that returns a vector of equations defining gate dynamics
end
```
"""
abstract type GateVarForm end 

struct AlphaBeta <: GateVarForm end 
struct SteadyStateTau <: GateVarForm end
struct SimpleGate <: GateVarForm end
struct ParameterGate <: GateVarForm end
struct HeavisideSum <: GateVarForm end

struct Gate{T<:GateVarForm} <: AbstractGatingVariable
    form::Type{T}
    output::Num
    eqs::Vector{Equation}
    props::Dict{Symbol,Any}
end

"""
$(TYPEDEF)

Low-level contructor for `Gate`.

A gate has a single symbolic `output` and stores properties (passed as a variable length
list of keyword arguments). Gate properties are accessible via `get` and `getproperty`.

# Example
```jldoctest; setup=:(using Conductor,ModelingToolkit; struct MyGateType<:Conductor.GateVarForm end)
julia> @variables t X(t)
2-element Vector{Num}:
    t
 X(t)

julia> g = Gate{MyGateType}(MyGateType, X, Equation[], prop1 = "foo", prop2 = 62)
Gate{MyGateType}(MyGateType, X(t), Symbolics.Equation[], Dict{Symbol, Any}(:prop2 => 62, :prop1 => "foo"))

julia> (g.prop1, g.prop2)
("foo", 62)
```
"""
function Gate{T}(form::Type{T}, output::Num, eqs::Vector{Equation}; kwargs...) where T <: GateVarForm
    return Gate{T}(form, output, eqs, kwargs)
end

# Internal API: Gate property getters
steadystate(x::AbstractGatingVariable) = get(x, :ss, nothing)
Base.exponent(x::AbstractGatingVariable) = get(x, :p, 1)
ModelingToolkit.get_eqs(x::AbstractGatingVariable, chan = nothing) = getfield(x, :eqs)

"""
$(TYPEDSIGNATURES)

Accepts expressions for forward (??) and reverse (??) reaction rates as descriptors for its
kinetics.

See also: [`get_eqs`](@ref).
"""
function Gate(form::Type{AlphaBeta}, ??, ??; name = Base.gensym("GateVar"), kwargs...)
    x??? = ??/(?? + ??)
    x = only(@variables $name(t) = x???)
    eqs = [D(x) ~ ??*(1 - x) - ??*x]
    return Gate{AlphaBeta}(form, x, eqs; ss = x???, kwargs...)
end

"""
$(TYPEDSIGNATURES)

Accepts expressions for its steady-state activation, x???(V???), and the time constant, ?????(V???),
as descriptors for its kinetics.

See also: [`get_eqs`](@ref).
"""
function Gate(form::Type{SteadyStateTau}, x???, ?????; name = Base.gensym("GateVar"), kwargs...)
    x = only(@variables $name(t) = x???)
    eqs = [D(x) ~ inv(?????)*(x??? - x)]
    return Gate{SteadyStateTau}(form, x, eqs; ss = x???, kwargs...)
end

function Base.convert(::Type{Gate{SimpleGate}},
                      x::Union{Gate{AlphaBeta},Gate{SteadyStateTau}})
    return Gate(SimpleGate, steadystate(x), p = exponent(x),
                name = Symbolics.tosymbol(output(x), escape=false))
end

"""
$(TYPEDSIGNATURES)

Accepts any symbolic expression as an explicit definition of the gate dynamics.
"""
function Gate(form::Type{SimpleGate}, rhs; default = rhs, name = Base.gensym("GateVar"), kwargs...)
    x = only(@variables $name(t) = default)
    return Gate{SimpleGate}(form, x, [x ~ rhs]; kwargs...)
end

"""
$(TYPEDSIGNATURES)

A static parameter gate with initial value, `val`.
"""
function Gate(form::Type{ParameterGate}, val; name = Base.gensym("GateVar"), kwargs...)
    x = only(@parameters $name = val)
    return Gate{ParameterGate}(form, x, Equation[]; val = val, kwargs...)
end

"""
    Gate(::Type{HeavisideSum}; threshold = 0mV, decay = 150, name = Base.gensym("GateVar") [, saturation])

Synaptically-activated dynamics. Sums the step-function values for presynaptic (extrinsic)
voltages.

The optional argument `saturation` sets a upper limit on the value of this gate.

See also: [`get_eqs`](@ref).
"""
function Gate(form::Type{HeavisideSum}; threshold = 0mV, decay = 150,
              name = Base.gensym("GateVar"), kwargs...) 
    x = only(@variables $name(t) = 0.0) # synaptically activated gate inits to 0.0
    return Gate{HeavisideSum}(form, x, Equation[]; threshold = threshold, decay = decay,
                              kwargs...)
end 

"""
    get_eqs(var::Gate{HeavisideSum}, chan)

Returns an equation of the form:

``
\\frac{dx}{dt}={\\displaystyle \\sum_{j}Heaviside(V_{j}-threshold)-x/saturation} 
``
"""
function ModelingToolkit.get_eqs(var::Gate{HeavisideSum}, chan)
    thold, decay = var.threshold, var.decay
    thold_val = ustrip(Float64, mV, thold)
    out = output(var)
    isempty(subscriptions(chan)) && return [D(out) ~ 0]
    @named V??? = ExtrinsicPotential(n = length(subscriptions(chan))) 
    # Derived from Pinsky & Rinzel 1994 - Equation 4 
    # S'??? = ??? ????(V??? - 10) - S???/150
    saturation = get(var, :saturation, nothing)
    if isnothing(saturation)
        return [D(out) ~ sum(V??? .>= thold_val) .- (out/decay)]
    else
        # out cannot continue to grow past the saturation limit
        return [D(out) ~ (out < saturation)*sum(V??? .>= thold_val) .- (out/decay)]
    end
end

#"""
#    get_eqs(var::Gate{<:Union{AlphaBeta, SteadySateTau}}, chan)
#
#Generate the voltage- and time-dependent differential equation modeling the output of a gate
#that was specified with [??(V???), ??(V???)] or [x???(V???), ?????(V???)].
#
#For `Gate{SteadyStateTau}`, the model form is:
#
#``
#\\frac{dx}{dt}=\\frac{1}{\\tau_{x}(x_{\\infty}-x)}
#``
#
#For `Gate{AlphaBeta}`, the model form is:
#
#``
#\\frac{dx}{dt} = \\alpha_{x}(1-x)-\\beta_{x}x
#``
#"""
#function ModelingToolkit.get_eqs(var::Gate{<:Union{AlphaBeta,SteadyStateTau}}, chan)
#    x, x???, ????? = output(var), steadystate(var), timeconstant(var)
#    return [D(x) ~ inv(?????)*(x??? - x)]
#end

#ModelingToolkit.get_eqs(var::Gate{SteadyState}, chan) = [output(var) ~ steadystate(var)]
#ModelingToolkit.get_eqs(var::Gate{ConstantValue}, chan) = Equation[]

############################################################################################
# Macros (needs updating)
############################################################################################
#=
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
    filter!(x -> x != :V???, param_list)

    if issetequal([:??, :??], rate_vars)  # Equation order doesn't matter
        rate_type = AlphaBeta
        gate_expr = :(Gate($AlphaBeta, ??, ??, $gate_pow; name=$(QuoteNode(name))))
    elseif issetequal([:ss, :??], rate_vars)  # Equation order doesn't matter
        rate_type = SteadyStateTau
        gate_expr = :(Gate($SteadyStateTau, ss, ??, $p.args[2]; name=$(QuoteNode(name))))
    else
        throw("invalid keywords")
    end
    pushfirst!(ex.args, :(V??? = MembranePotential()))
    foreach(param -> pushfirst!(ex.args, :(@parameters $param)), param_list)  # generate parameters
    push!(ex.args, gate_expr)
    ex
end
=#
