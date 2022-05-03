
using Graphs

abstract type AbstractSynapse{T} <: AbstractEdge{T} end
abstract type AbstractNeuralCircuitGraph{T} <: AbstractGraph{T} end

struct Synapse{T} <: AbstractSynapse{T}
    source::CompartmentSystem
    target::CompartmentSystem
    syntype::ConductanceSystem
end

pre(x::Synapse) = getfield(x, :source)
post(x::Synapse) = getfield(x, :target)
class(x::Synapse) = getfield(x, :syntype)

# Synapse(neuron1, neuron2, Glut, 100pS)
# Expect a list of eltype Synapse => construct topology as Dict{SynapticChannel, MetaDiGraph}

abstract type AbstractNeuralCircuitSystem <: AbstractTimeDependentSystem end

struct NeuralCircuit <: AbstractNeuralCircuitSystem
    ivs::Num
    topology::Dict # implement a multiplexed network using a dictionary to access layers
    extensions::Vector{ODESystem}
    defaults::Dict
    name::Symbol
end

get_topology(x::AbstractNeuralCircuitSystem) = getfield(x, :topology)

function NeuralCircuit(graph, extensions::Vector{ODESystem} = [];
                       defaults = Dict(), name::Symbol = Base.gensym(:Network))
    
    neurons = Set()
    synapse_types = Set()
    topology = Dict()

    for synapse in synapses
        push!(neurons, pre(synapse))
        push!(neurons, post(synapse))
        push!(synapse_types, nameof(class(synapse)))
    end
    
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
