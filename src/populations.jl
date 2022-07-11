
struct Population{T<:AbstractCompartmentSystem}
    "The number of neurons in the population." 
    N::Int
    "The base neuron model."
    proto::T
    "Name of the population."
    popname::Symbol
end

Base.nameof(P::Population) = getfield(P, :popname)
Base.eltype(::Population{T}) where {T} = T
Base.length(P::Population) = getfield(P, :N) 

function Base.iterate(P::Population, state=1)
    state > P.N && return nothing
    neuron = MTK.rename(P.proto, Symbol(nameof(P), '₊', nameof(P.proto), '_', state))
    return (neuron, state+1)
end

function Base.iterate(rP::Iterators.Reverse{Population{T}}, state=rP.itr.N) where {T<:AbstractCompartmentSystem}
    state < 1 && return nothing
    pop = rP.itr
    neuron = MTK.rename(pop.proto, Symbol(nameof(pop), '₊', nameof(pop.proto), '_', state))
    return (neuron, state-1)
end

function Base.getindex(P::Population, i::Int)
    1 <= i <= P.N || throw(BoundsError(P, i))
    return MTK.rename(P.proto, Symbol(nameof(P), '₊', nameof(P.proto), '_', i))
end

Base.firstindex(P::Population) = 1
Base.lastindex(P::Population) = length(P)
Base.getindex(P::Population, i::Number) = P[convert(Int, i)]
Base.getindex(P::Population, I) = [P[i] for i in I]


