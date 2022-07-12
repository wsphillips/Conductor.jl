"""
$(TYPEDEF)

An iterable population (or small group) of neurons.

A `Population` is a lazy representation of neurons templated from a single model prototype.
Iterating over a `Population` returns up to _n_ uniquely-named copies of the prototype that
are serially numbered and, optionally, namespaced. Parameter values of the replicate neurons
can be individually defined by indexed assignment, or applied to the entire group via
distributions.

$(FIELDS)
"""
struct Population{T<:AbstractCompartmentSystem}
    "The base neuron model."
    proto::T
    "The number of neurons in the population." 
    n::Int
    "Name of the population."
    popname::Symbol
    "Population defaults"
    defs::Vector{Dict{Num, Any}}
    "Stimuli applied to select neurons."
    stimuli::Vector{Vector{Equation}}
end

function Population(neuron, n; name = Symbol(""))
    Population(neuron, n, name, [Dict{Num, Any}() for _ in 1:n], [Equation[] for _ in 1:n])
end

Base.nameof(p::Population) = getfield(p, :popname)
Base.eltype(::Population{T}) where {T} = T
Base.length(p::Population) = getfield(p, :n) 
prototype(p::Population) = getfield(p, :proto)

function Base.iterate(p::Population{T}, state=1) where {T}
    state > length(p) && return nothing
    new_name = Symbol(nameof(p), '₊', nameof(prototype(p)), '_', state)
    neuron = T(prototype(p); name = new_name, defaults = p.defs[state])
    return (neuron, state+1)
end

function Base.iterate(rP::Iterators.Reverse{Population{T}}, state=length(rP.itr)) where {T<:AbstractCompartmentSystem}
    state < 1 && return nothing
    pop = rP.itr
    new_name = Symbol(nameof(pop), '₊', nameof(prototype(pop)), '_', state)
    neuron = T(prototype(p); name = new_name, defaults = p.defs[state])
    return (neuron, state-1)
end

function Base.getindex(p::Population, i::Int)
    1 <= i <= length(p) || throw(BoundsError(p, i))
    return MTK.rename(prototype(p), Symbol(nameof(p), '₊', nameof(prototype(p)), '_', i))
end

Base.firstindex(p::Population) = 1
Base.lastindex(p::Population) = length(p)
Base.getindex(p::Population, i::Number) = p[convert(Int, i)]
Base.getindex(p::Population, I) = [p[i] for i in I]

function (pop::Population)(p::Pair{Num,D}) where {D<:Distribution}
    new_defs = [Dict{Num, Any}() for _ in 1:length(pop)]
    for neuron_defs in new_defs
        push!(neuron_defs, p.first => rand(p.second))
    end
    return @set pop.defs = merge.(pop.defs, new_defs)
end

function (pop::Population)(p::Pair{Num,Q}) where {Q<:Quantity}
    var, quant = p.first, p.second
    cunits = hasmetadata(var, ConductorUnits) ? getmetadata(var, ConductorUnits) : nothing
    isnothing(cunits) && throw("Units do not match.")
    val = ustrip(Float64, cunits, quant)
    new_defs = [Dict{Num, Any}() for _ in 1:length(pop)]
    for neuron_defs in new_defs
        push!(neuron_defs, var => val)
    end
    return @set pop.defs = merge.(pop.defs, new_defs)
end

function Base.setindex!(pop::Population, p::Pair{Num,D}, i::Int) where {D<:Distribution}
    var, dist = p.first, p.second
    push!(pop.defs[i], var => rand(dist))
end

function Base.setindex!(pop::Population, p::Pair{Num,Q}, i::Int) where {Q<:Quantity}
    var, quant = p.first, p.second
    cunits = hasmetadata(var, ConductorUnits) ? getmetadata(var, ConductorUnits) : nothing
    isnothing(cunits) && throw("Units do not match.")
    val = ustrip(Float64, cunits, quant)
    push!(pop.defs[i], var => val)   
end

