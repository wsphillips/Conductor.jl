
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
    function NeuronalNetworkSystem(iv, synapses, extensions, name, eqs, systems; checks = false)
        if checks
            # placeholder
        end
        new(iv, synapses, extensions, name, eqs, systems)
    end
end

function NeuronalNetworkSystem(synapses::Vector{Synapse}, extensions::Vector{AbstractTimeDependentSystem} = AbstractTimeDependentSystem[];
                               name::Symbol = Base.gensym(:Network))
    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    return NeuronalNetworkSystem(t, synapses, extensions, name, eqs, systems)
end

#get_topology(x::AbstractNeuronalNetworkSystem) = getfield(x, :topology)
get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)
get_synapses(x::AbstractNeuronalNetworkSystem) = getfield(x, :synapses)

# Forward getters to internal system
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
    union!(getfield(x, :systems), build_toplevel(x)[5])
end

function Base.convert(::Type{ODESystem}, x::NeuronalNetworkSystem) end

function build_toplevel!(dvs, ps, eqs, defs, network_sys::NeuronalNetworkSystem)
    
    synapses = get_synapses(network_sys)
    preneurons = Set()
    postneurons = Set()
    synapse_types = Set()
    reversals = Set()
    
    voltage_fwds = Set{Equation}()
    
    # Bin the fields
    for synapse in synapses
        push!(preneurons, presynaptic(synapse))
        push!(preneurons, postsynaptic(synapse))
        push!(synapse_types, class(synapse))
        push!(reversals, reversal(synapse))
    end
    
    # Reset all synaptic information
    allneurons = union(preneurons, postneurons)
    foreach(x -> empty!(get_synapses(x)), allneurons)
    foreach(x -> empty!(get_synaptic_reversals(x)), allneurons)

    # Push synaptic information to each neuron
    for synapse in synapses
        syn = replicate(class(synapse))
        post = postsynaptic(synapse)
        pre  = presynaptic(synapse)
        push!(get_synapses(post), syn)
        push!(get_synaptic_reversals(post), reversal(synapse))
        # Note this will probably cause an error...
        push!(voltage_fwds, getproperty(post, nameof(syn)).Vₘ ~ pre.Vₘ) 
    end

    union!(eqs, voltage_fwds)
    return dvs, ps, eqs, defs, collect(allneurons)
end

#=
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
