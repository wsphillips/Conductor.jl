
"""
AbstractSynapse

    -must have at least fields: source, target, conductance
    - eventually subtype AbstractEdge and implement Graphs interface?
"""
abstract type AbstractSynapse end # <: AbstractEdge{T} end

struct Synapse <: AbstractSynapse
    source::CompartmentSystem
    target::CompartmentSystem
    conductance::AbstractConductanceSystem
    reversal::Num
end

function Synapse(x::Pair, cond, rev)
    return Synapse(x.first, x.second, cond, rev) 
end

# Utility function to copy a system, but generate a new name (gensym)
# Perhaps we could cache a global index for synapse types? Then we could
# rename using a plain integer for better readability.
function replicate(x::Union{AbstractCompartmentSystem,AbstractConductanceSystem})
    rootname = ModelingToolkit.getname(x)
    new = deepcopy(x)
    return ModelingToolkit.rename(new, Base.gensym(rootname))
end

presynaptic(x::Synapse) = getfield(x, :source)
postsynaptic(x::Synapse) = getfield(x, :target)
class(x::Synapse) = getfield(x, :conductance)
reversal(x::Synapse) = getfield(x, :reversal)
abstract type AbstractNeuronalNetworkSystem <: AbstractTimeDependentSystem end

struct NeuronalNetworkSystem <: AbstractNeuronalNetworkSystem
    #topology::AbstractNetworkTopology
    iv::Num
    synapses::Vector{AbstractSynapse}
    extensions::Vector{AbstractTimeDependentSystem}
    name::Symbol
    eqs::Vector{Equation}
    systems::Vector{AbstractTimeDependentSystem}
    observed::Vector{Equation}
    function NeuronalNetworkSystem(iv, synapses, extensions, name, eqs, systems, observed; checks = false)
        if checks
            # placeholder
        end
        new(iv, synapses, extensions, name, eqs, systems, observed)
    end
end

function NeuronalNetworkSystem(synapses::Vector{Synapse}, extensions::Vector{AbstractTimeDependentSystem} = AbstractTimeDependentSystem[];
                               name::Symbol = Base.gensym(:Network))
    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    return NeuronalNetworkSystem(t, synapses, extensions, name, eqs, systems, observed)
end

#get_topology(x::AbstractNeuronalNetworkSystem) = getfield(x, :topology)
get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)
get_synapses(x::AbstractNeuronalNetworkSystem) = getfield(x, :synapses)

MTK.get_states(x::AbstractNeuronalNetworkSystem) = collect(build_toplevel(x)[1])
MTK.has_ps(x::AbstractNeuronalNetworkSystem) = !isempty(build_toplevel(x)[2])
MTK.get_ps(x::AbstractNeuronalNetworkSystem) = collect(build_toplevel(x)[2])

function MTK.get_eqs(x::AbstractNeuronalNetworkSystem)
    empty!(getfield(x, :eqs))
    union!(getfield(x, :eqs), build_toplevel(x)[3])
end

MTK.get_defaults(x::AbstractNeuronalNetworkSystem) = build_toplevel(x)[4]

function MTK.get_systems(x::AbstractNeuronalNetworkSystem)
    empty!(getfield(x, :systems))
    union!(getfield(x, :systems), build_toplevel(x)[5], get_extensions(x))
end

function Base.convert(::Type{ODESystem}, nnsys::NeuronalNetworkSystem)
    states, params, eqs, defs, allneurons = build_toplevel(nnsys)
    all_systems = map(x -> convert(ODESystem, x), allneurons)
    odesys = ODESystem(eqs, t, states, params; defaults = defs, name = nameof(nnsys))
    return compose(odesys, all_systems)
end

function build_toplevel!(dvs, ps, eqs, defs, network_sys::NeuronalNetworkSystem)
    
    synapses = get_synapses(network_sys)
    preneurons  = Set()
    postneurons = Set()
    neurons = Set()
    
    voltage_fwds = Set{Equation}()
    # Bin the fields
    for synapse in synapses
        push!(preneurons, presynaptic(synapse))
        push!(postneurons, postsynaptic(synapse))
    end
    
    union!(neurons, preneurons, postneurons)

    # Reset all synaptic information
    foreach(x -> empty!(get_synapses(x)), neurons)
    foreach(x -> empty!(get_synaptic_reversals(x)), neurons)

    # Push synaptic information to each neuron
    for synapse in synapses
        syn = replicate(class(synapse))
        #syn = class(synapse)
        post = postsynaptic(synapse)
        pre  = presynaptic(synapse)
        push!(get_synapses(post), syn)
        push!(get_synaptic_reversals(post), reversal(synapse))
        # Note this will probably cause an error...
        push!(voltage_fwds, pre.Vₘ ~ getproperty(post, nameof(syn)).Vₘ) 
    end

    union!(eqs, voltage_fwds)
    return dvs, ps, eqs, defs, collect(neurons)
end

#=
# Early scratch for graph based objects/network construction
using Graphs
using SparseArrays

function NeuronalNetworkSystem(topology::NetworkTopology, extensions::Vector{ODESystem} = [];
                         defaults = Dict(), name::Symbol = Base.gensym(:Network))
    
   
    for layer in synapse_types
        topology[layer] = MetaGraph(SimpleDiGraph(),
                                    VertexMeta = AbstractCompartmentSystem,
                                    EdgeMeta = AbstractSynapse) 
                         # weightfunction can be set for default weight value getter
        # vertices are always the same
        for neuron in neurons
            topology[layer][nameof(neuron)] = deepcopy(neuron)
        end
    end
    
    for synapse in synapses
        topology[syntype][nameof(pre(synapse)), nameof(post(synapse))] = deepcopy(synapse)
    end
    
    NeuralCircuit(t, topology, extensions, defaults, name)
end

"""
AbstractNeuronalNetworkTopology

    - multilayer graph representation of the network
    - eventually subtype AbstractGraph and implement Graphs interface?
"""
abstract type AbstractNetworkTopology end #{T} <: AbstractGraph{T} end

# note the eltype T _could_ be Num if we need it to
struct NetworkTopology{V<:AbstractCompartmentSystem,T} <: AbstractNetworkTopology
    layered_adjacency_matrix::Vector{SparseMatrixCSC{T,Int64}} # vector idx == layer
    vert_idx::IdDict{V,T} # lookup neuron index
    flayer_idx::IdDict{V,T} # forward lookup index from synapse type
    blayer_idx::IdDict{V,T} # backward lookup synapse type from index
end

function NetworkToplogy(synapses) end

function Base.getindex(x::NetworkTopology, pre::AbstractCompartmentSystem, post::AbstractCompartmentSystem)
    #getkey(collection, 
    return x.layered_adjacency_matrix[] 
end
=#
