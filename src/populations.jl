"""
$(TYPEDEF)

A population (or small group) of neurons.

A `Population` is a lazy representation of neurons templated from a single model prototype.
Iterating over a `Population` returns up to _n_ uniquely-named copies of the prototype that
are serially numbered and, optionally, namespaced. Parameter values of the replicate neurons
can be individually defined by indexed assignment, or applied to the entire group via
distributions.

$(FIELDS)
"""
struct Population{T <: AbstractCompartmentSystem}
    neurons::Vector{T}
    "Base name of the population."
    popname::Symbol
    "Population defaults"
    defs::Vector{Dict{Num, Any}}
    "Stimuli applied to select neurons."
    stimuli::Vector{Vector{Equation}}
end

# Until MTK supports equivalent systems (loops), we scalarize the population on construction.
function Population(neuron, n; name)
    neurons = map(i -> rename(deepcopy(neuron), Symbol(name, :_, i)), 1:n) #same convention as @named
    Population(neurons, name, [Dict{Num, Any}() for _ in 1:n], [Equation[] for _ in 1:n])
end

Base.nameof(p::Population) = getfield(p, :popname)
Base.eltype(::Population{T}) where {T} = T
Base.length(p::Population) = length(getfield(p, :neurons))

function Base.iterate(p::Population, state = 1)
    state > length(p) && return nothing
    neuron = p.neurons[state]
    return (neuron, state + 1)
end

function Base.iterate(rP::Iterators.Reverse{Population},
                      state = length(rP.itr))
    state < 1 && return nothing
    pop = rP.itr
    neuron = pop.neurons[state]
    return (neuron, state - 1)
end

function Base.getindex(p::Population, i::Int)
    1 <= i <= length(p) || throw(BoundsError(p, i))
    return p.neurons[i]
end

Base.firstindex(p::Population) = 1
Base.lastindex(p::Population) = length(p)
Base.getindex(p::Population, i::Number) = p[convert(Int, i)]
Base.getindex(p::Population, I) = [p[i] for i in I]

############################################################################################

# Setting specific distribution values with pair syntax
function (pop::Population)(p::Pair{Num, D}) where {D <: Distribution}
    new_defs = [Dict{Num, Any}() for _ in 1:length(pop)]
    for neuron_defs in new_defs
        push!(neuron_defs, p.first => rand(p.second))
    end
    return @set pop.defs = merge.(pop.defs, new_defs)
end

# Setting specific default values with pair syntax
function (pop::Population)(p::Pair{Num, Q}) where {Q <: Quantity}
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

# set _distribution_ of a symbol in the population (?)
function Base.setindex!(pop::Population, p::Pair{Num, D}, i::Int) where {D <: Distribution}
    var, dist = p.first, p.second
    push!(pop.defs[i], var => rand(dist))
end

# set explicit default of a symbol in the population
function Base.setindex!(pop::Population, p::Pair{Num, Q}, i::Int) where {Q <: Quantity}
    var, quant = p.first, p.second
    cunits = hasmetadata(var, ConductorUnits) ? getmetadata(var, ConductorUnits) : nothing
    isnothing(cunits) && throw("Units do not match.")
    val = ustrip(Float64, cunits, quant)
    push!(pop.defs[i], var => val)
end
