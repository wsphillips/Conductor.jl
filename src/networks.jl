
#using Graphs
    using SparseArrays
"""
AbstractSynapse

    -must have at least fields: source, target, conductance
    - eventually subtype AbstractEdge and implement Graphs interface?
"""
abstract type AbstractSynapse end # <: AbstractEdge{T} end

struct Synapse{T<:AbstractConductanceSystem} <: AbstractSynapse
    source::CompartmentSystem
    target::CompartmentSystem
    conductance::T
    reversal::Num
end

function Synapse(x::Pair, cond, rev)
    return Synapse(x.first, x.second, cond, rev) 
end

# Utility function to copy a system, but generate a new name (gensym)
function replicate(x::Union{AbstractCompartmentSystem,AbstractConductanceSystem})
    rootname = ModelingToolkit.getname(x)
    new = deepcopy(x)
    return ModelingToolkit.rename(new, Base.gensym(rootname))
end

presynaptic(x::Synapse) = getfield(x, :source)
postsynaptic(x::Synapse) = getfield(x, :target)
class(x::Synapse) = getfield(x, :conductance)
equilibrium(x::Synapse) = getfield(x, :reversal)
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

function NetworkToplogy(synapses)

    neurons = Set()
    synapse_classes = Set()

    for synapse in synapses
        push!(neurons, presynaptic(synapse))
        push!(neurons, postsynaptic(synapse))
        # using nameof is a bad heuristic
        push!(synapse_classes, nameof(class(synapse)))
    end

    return
end

#=
function Base.getindex(x::NetworkTopology, pre::AbstractCompartmentSystem, post::AbstractCompartmentSystem)
    # getkey(collection, 
    return x.layered_adjacency_matrix[] 
end
=#

# Synapse(neuron1, neuron2, Glut, 100pS)
# Expect a list of eltype Synapse => construct topology as Dict{SynapticChannel, MetaDiGraph}

abstract type AbstractNeuronalNetworkSystem <: AbstractTimeDependentSystem end

struct NeuronalNetworkSystem{S<:AbstractTimeDependentSystem} <: AbstractNeuronalNetworkSystem
    topology::AbstractNetworkTopology
    extensions::Vector{ODESystem}
    sys::S
    name::Symbol
end

get_topology(x::AbstractNeuronalNetworkSystem) = getfield(x, :topology)
get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)

Base.convert(::Type{ODESystem}, x::NeuronalNetworkSystem{ODESystem}) = getfield(x, :sys)

# Forward getters to internal system
MTK.get_systems(x::AbstractNeuronalNetworkSystem) = get_systems(getfield(x, :sys))
MTK.get_eqs(x::AbstractNeuronalNetworkSystem) = get_eqs(getfield(x, :sys))
MTK.get_dvs(x::AbstractNeuronalNetworkSystem) = get_dvs(getfield(x, :sys))
MTK.has_ps(x::AbstractNeuronalNetworkSystem) = MTK.has_ps(getfield(x, :sys))
MTK.get_ps(x::AbstractNeuronalNetworkSystem) = get_ps(getfield(x, :sys))
MTK.get_defaults(x::AbstractNeuronalNetworkSystem) = get_defaults(getfield(x, :sys))
MTK.get_states(x::AbstractNeuronalNetworkSystem) = get_states(getfield(x, :sys))
MTK.get_ivs(x::AbstractNeuronalNetworkSystem) = get_ivs(getfield(x, :sys))
MTK.get_iv(x::AbstractNeuronalNetworkSystem) = get_iv(getfield(x, :sys))
MTK.independent_variables(x::AbstractNeuronalNetworkSystem) = MTK.independent_variables(getfield(x, :sys))
MTK.get_observed(x::AbstractNeuronalNetworkSystem) = MTK.get_observed(getfield(x, :sys))

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

function NeuronalNetworkSystem(synapses::Vector{Synapse}, extensions::Vector{ODESystem} = [];
                         defaults = Dict(), name::Symbol = Base.gensym(:Network))

    neurons = Set()
    synapse_types = Set()

    for synapse in synapses
        push!(neurons, pre(synapse))
        push!(neurons, post(synapse))
        push!(synapse_types, nameof(class(synapse)))
    end

    return NeuronalNetwork(topology, extensions, defaults = defaults, name = name)
end

#=
function get_eqs(x::AbstractNeuralCircuitSystem) end
function get_states(x::AbstractNeuralCircuitSystem) end
MTK.has_ps(x::NeuralCircuit) = any(MTK.has_ps, getfield(x, :extensions))
function get_ps(x::AbstractNeuralCircuitSystem) end
function defaults(x::AbstractNeuralCircuitSystem) end
function get_systems(x::AbstractNeuralCircuitSystem) end
function Base.convert(::Type{ODESystem}, network::NeuralCircuit) end

function get_synapses(x::NeuralCircuit, layer::Symbol)
    keys(getfield(get_topology(x)[layer], :eprops))
end

function get_neurons(x::NeuralCircuit)
    top = get_topology(x)
    g = first(first(top)) # grab the first layer
    neurons = Set(values()) # FINISH ME
end

function build_toplevel!(dvs, ps, eqs, defs, x::NeuralCircuit)
    
    # Clear existing synaptic conductances and synaptic reversals from neurons
    # use empty!(set)
end

function build_toplevel(x::NeuralCircuit)
    dvs
    ps
    eqs
    defs
    build_toplevel!(dvs, ps, eqs, defs, x)
    return dvs, ps, dqs, defs
end

## Old constructor below

function NeuralCircuit(neurons, topology; name = Base.gensym(:Network))

    all_neurons = Set(getproperty.(neurons, :sys))
    eqs = Equation[]
    params = Set()
    states = Set() # place holder
    defaultmap = Pair[]
    all_systems = AbstractSystem[]
    push!(all_systems, all_neurons...)
    post_neurons = Set()

    # what types of synapses do we have
    all_synapses = [synapse.second[2] for synapse in topology]
    all_synapse_types = [nameof(x.sys) for x in all_synapses]
    synapse_types = unique(all_synapse_types)

    # how many of each synapse type are there
    synapse_counts = Dict{Symbol,Int64}()

    for n in synapse_types
        c = count(x -> isequal(x, n), all_synapse_types)
        push!(synapse_counts, n => c)
    end

    # Extract reversal potentials
    reversals = unique([x.reversal for x in all_synapses])
    push!(params, reversals...)
    rev_meta = [getmetadata(x, ConductorEquilibriumCtx).val for x in reversals]

    for (i,j) in zip(reversals, rev_meta)
        push!(defaultmap, i => ustrip(Float64, mV, j))
    end

    voltage_fwds = Set()

    # create a unique synaptic current for each post-synaptic target
    for synapse in topology
        pre = synapse.first.sys # pre-synaptic neuron system
        post = synapse.second[1].sys # post-synaptic neuron system
        push!(post_neurons, post)
        Erev = synapse.second[2].reversal
        syntype = synapse.second[2].sys # synapse system
        synname = nameof(syntype)
        syn = @set syntype.name = Symbol(synname, synapse_counts[synname]) # each synapse is uniquely named
        synapse_counts[synname] -= 1
        push!(all_systems, syn)
        push!(voltage_fwds, syn.Vₘ ~ pre.Vₘ)

        if post.Isyn ∈ Set(vcat((get_variables(x.lhs) for x in eqs)...))
            idx = findfirst(x -> isequal([post.Isyn], get_variables(x.lhs)), eqs)
            eq = popat!(eqs, idx)
            eq = eq.lhs ~ eq.rhs + (syn.g * (post.Vₘ - Erev))
            push!(eqs, eq)
        else
            push!(eqs, post.Isyn ~ syn.g * (post.Vₘ - Erev))
        end
    end

    for nonpost in setdiff(all_neurons, post_neurons)
        push!(eqs, D(nonpost.Isyn) ~ 0)
    end

    append!(eqs, collect(voltage_fwds))
    network_system = ODESystem(eqs, t, states, params; defaults = defaultmap, name = name )
    return compose(network_system, all_systems)
end
=#
