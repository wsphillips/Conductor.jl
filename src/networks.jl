
# gets converted into an edge type for MetaGraphs
struct Synapse
    source
    target
    conductance_type
    weight
end

pre(x::Synapse) = getfield(x, :source)
post(x::Synapse) = getfield(x, :target)
class(x::Synapse) = getfield(x, :conductance_type)
weight(x::Synapse) = getfield(x, :weight)

# Synapse(neuron1, neuron2, Glut, 100pS)
# Expect a list of eltype Synapse => construct topology as Dict{SynapticChannel, MetaDiGraph}

abstract type AbstractNetworkSystem <: AbstractTimeDependentSystem

struct NetworkSystem <: AbstractNetworkSystem end
    ivs::Num
    topology::Dict # implement a multiplexed network using a dictionary to access layers
    extensions::Vector{ODESystem}
    defaults::Dict
    name::Symbol
end

function NetworkSystem(synapses, extensions::Vector{ODESystem} = []; defaults = Dict(), name::Symbol = Base.gensym(:Network))
    
    neurons = Set()
    synapse_types = Set()
    topology = Dict()

    for synapse in synapses
        push!(neurons, pre(synapse))
        push!(neurons, post(synapse))
        push!(synapse_types, class(synapse))
    end
    
    for layer in synapse_types
        topology[layer] = MetaGraph(SimpleDiGraph(),
                                      VertexMeta = AbstractCompartmentSystem,
                                      EdgeMeta = AbstractConductanceSystem) 
                         # weightfunction can be set for default weight value getter
        for neuron in neurons
            topology[layer][nameof(neuron)] = neuron
        end
    end
    
    for synapse in synapses
        syntype = class(synapse)
        wt = weight(synapse)
        topology[class(synapse)][nameof(pre(synapse)), nameof(post(synapse))] = syntype(wt)
    end
    
    NetworkSystem(t, topology, extensions, defaults, name)
end

function build_toplevel!(dvs, ps, eqs, defs, x::NetworkSystem) end

function build_toplevel(x::NetworkSystem)
    dvs
    ps
    eqs
    defs
    build_toplevel!(dvs, ps, eqs, defs, x)
    return dvs, ps, dqs, defs
end

function get_eqs(x::AbstractNetworkSystem) end

function get_states(x::AbstractNetworkSystem) end

MTK.has_ps(x::NetworkSystem) = any(MTK.has_ps, getfield(x, :extensions))

function get_ps(x::AbstractNetworkSystem) end

function defaults(x::AbstractNetworkSystem) end

function get_systems(x::AbstractNetworkSystem) end

function Base.convert(::Type{ODESystem}, network::NetworkSystem) end

function NetworkSystem(neurons, topology; name = Base.gensym(:Network))

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
