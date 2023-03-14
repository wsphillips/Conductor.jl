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

struct Gate{T <: GateVarForm} <: AbstractGatingVariable
    form::Type{T}
    output::Num
    eqs::Vector{Equation}
    props::Dict{Symbol, Any}
end

"""
$(TYPEDEF)

Low-level constructor for `Gate`.

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
function Gate{T}(form::Type{T}, output::Num, eqs::Vector{Equation};
                 kwargs...) where {T <: GateVarForm}
    return Gate{T}(form, output, eqs, kwargs)
end

# Internal API: Gate property getters
steadystate(x::AbstractGatingVariable) = get(x, :ss, nothing)
Base.exponent(x::AbstractGatingVariable) = get(x, :p, 1)
ModelingToolkit.get_eqs(x::AbstractGatingVariable, chan = nothing) = getfield(x, :eqs)

"""
$(TYPEDSIGNATURES)

Accepts expressions for forward (α) and reverse (β) reaction rates as descriptors for its
kinetics.

See also: `get_eqs`.
"""
function Gate(form::Type{AlphaBeta}, α, β; name = Base.gensym("GateVar"), kwargs...)
    x∞ = α / (α + β)
    x = only(@variables $name(t)=x∞ [unit = NoUnits])
    eqs = [D(x) ~ (α * (1 - x) - β * x)]
    return Gate{AlphaBeta}(form, x, eqs; ss = x∞, kwargs...)
end

"""
$(TYPEDSIGNATURES)

Accepts expressions for its steady-state activation, x∞(Vₘ), and the time constant, τₓ(Vₘ),
as descriptors for its kinetics.

See also: `get_eqs`.
"""
function Gate(form::Type{SteadyStateTau}, x∞, τₓ; name = Base.gensym("GateVar"), kwargs...)
    x = only(@variables $name(t)=x∞ [unit = NoUnits])
    eqs = [D(x) ~ inv(τₓ) * (x∞ - x)]
    return Gate{SteadyStateTau}(form, x, eqs; ss = x∞, kwargs...)
end

function Base.convert(::Type{Gate{SimpleGate}},
                      x::Union{Gate{AlphaBeta}, Gate{SteadyStateTau}})
    return Gate(SimpleGate, steadystate(x), p = exponent(x),
                name = Symbolics.tosymbol(output(x), escape = false))
end

"""
$(TYPEDSIGNATURES)

Accepts any symbolic expression as an explicit definition of the gate dynamics.
"""
function Gate(form::Type{SimpleGate}, rhs::Num; default = rhs, name = Base.gensym("GateVar"),
              kwargs...)
    x = only(@variables $name(t)=default [unit = NoUnits])
    return Gate{SimpleGate}(form, x, [x ~ rhs]; kwargs...) # rhs must return NoUnits
end

function Gate(form::Type{SimpleGate}, output::Num, eqs, kwargs...)
    Gate{SimpleGate}(form, output, eqs, Dict(kwargs))
end

Gate(rhs::Num; default = rhs, name = Base.gensym("GateVar"), kwargs...) = 
Gate(SimpleGate, rhs; default, name, kwargs...)

Gate(output::Num, eqs::Vector{Equation}, kwargs...) = Gate(SimpleGate, output, eqs, kwargs...)


"""
$(TYPEDSIGNATURES)

A static parameter gate with initial value, `val`.
"""
function Gate(form::Type{ParameterGate}, val; name = Base.gensym("GateVar"), kwargs...)
    x = only(@parameters $name=val [unit = NoUnits])
    return Gate{ParameterGate}(form, x, Equation[]; val = val, kwargs...)
end

