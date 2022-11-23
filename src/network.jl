
struct NetworkTopology
    multigraph::Dict{Any, SparseMatrixCSC{Float64, Int64}}
    neurons::Vector{Vector{CompartmentSystem}}
end

synaptic_models(topology::NetworkTopology) = keys(graph(topology))

function NetworkTopology(neurons, synaptic_models)
    m = length(neurons)
    n = sum(length(neuron) for neuron in neurons)
    multigraph = Dict(x => sparse(Int64[], Int64[], Float64[], m, n) for x in synaptic_models)
    return NetworkTopology(multigraph, neurons)
end

function NetworkTopology(g::SimpleDiGraph, neurons, synaptic_model,
                         default_weight = get_gbar(synaptic_model))
    multigraph = Dict(synaptic_model => adjacency_matrix(g) * default_weight)
    return NetworkTopology(multigraph, neurons)
end

neurons(topology::NetworkTopology) = getfield(topology, :neurons)
compartments(topology::NetworkTopology) = vcat(neurons(topology))
root_compartments(topology::NetworkTopology) = [first(neuron) for neuron in neurons(topology)]
graph(topology::NetworkTopology) = getfield(topology, :multigraph)

function add_synapse!(topology, pre, post, synaptic_model, weight)
    src = find_source(pre, topology)
    dst = find_target(post, topology)
    g = graph(topology)[synaptic_model]
    g[src, dst] = weight
end

function remove_synapse!(topology, pre, post, synaptic_model)
    src = find_source(pre, topology)
    dst = find_target(post, topology)
    g = graph(topology)[synaptic_model]
    g[src, dst] = zero(Num)
    dropzeros!(g)
end

# FIXME: needs update
function add_layer!(topology, synaptic_model::ConductanceSystem,
                    g = SimpleDiGraph(length(vertices(topology))),
                    default_weight = get_gbar(synaptic_model))
    push!(graph(topology), synaptic_model => adjacency_matrix(g) * default_weight)
end

Base.eltype(nt::NetworkTopology) = AbstractCompartmentSystem
Base.length(nt::NetworkTopology) = length(neurons(nt))

function Base.iterate(nt::NetworkTopology, state = 1)
    state > length(nt) && return nothing
    return (neurons(nt)[state], state + 1)
end

function Base.iterate(rev_nt::Iterators.Reverse{NetworkTopology},
                      state = length(rev_nt.itr))
    state < 1 && return nothing
    nt = rev_nt.itr
    return (neurons(nt)[state], state - 1)
end

function Base.getindex(nt::NetworkTopology, i)
    return neurons(nt)[i]
end

Base.firstindex(nt::NetworkTopology) = 1
Base.lastindex(nt::NetworkTopology) = length(nt)

function Base.setindex!(nt::NetworkTopology, g::SimpleDiGraph, cond::ConductanceSystem)
    add_layer!(nt, cond, g)
end

function Base.setindex!(nt::NetworkTopology, cond::ConductanceSystem, pre::T1,
                        post::T2) where {T1 <: AbstractCompartmentSystem,
                                         T2 <: AbstractCompartmentSystem}
    add_synapse!(nt, pre, post, cond)
end

function Base.setindex!(nt::NetworkTopology, weight::Real, pre, post)
    add_synapse!(nt, pre, post, weight)
end

"""
$(TYPEDEF)

A network of neurons with synaptic connections.

$(TYPEDFIELDS)
"""
struct NeuronalNetworkSystem <: AbstractNeuronalNetworkSystem
    # MTK fields
    eqs::Vector{Equation}
    "Independent variabe. Defaults to time, ``t``."
    iv::Num
    states::Vector
    ps::Vector
    observed::Vector{Equation}
    name::Symbol
    systems::Vector{AbstractTimeDependentSystem}
    defaults::Dict
    # Conductor fields
    topology::NetworkTopology
    reversal_map::Dict
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    function NeuronalNetworkSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                                   topology,
                                   reversal_map,
                                   extensions; checks = false)
        if checks
            # placeholder
        end
        new(eqs, iv, states, ps, observed, name, systems, defaults, topology, reversal_map,
            extensions)
    end
end

synaptic_models(sys::NeuronalNetworkSystem) = synaptic_models(topology(sys))

# events are handled by callbacks and do not modify
function connect_synapses!(gen, syn_model, compartments, topology)
    return compartments
end

function connect_synapses!(gen, syn_model::Union{CurrentSystem{T},ConductanceSystem{T}},
                           compartments, topology) where {T<:IntegratedSynapse}

    (; eqs, systems, defs) = gen
    new_compartments = deepcopy(compartments)
    roots = root_compartments
    g = multigraph[syn_model]   
    rows = rowvals(g)
    vals = nonzeros(g)
    voltage_fwds = Set{Equation}()
    
    # Synaptic conductance models are duplicated here and individually connected, relying
    # on users to specify extrinsic voltage symbols. In the future, when symbolic arrays and
    # component types become available we will probably rewrite how this works.
    for i in axes(g, 2) # for each set of presynaptic neurons per postsynaptic neuron
        model_copies = []
        pre_indexes = nzrange(g, i)
        isempty(pre_indexes) && continue # skip neurons with no synaptic input

        post_compartment = compartments[i]
        pre_compartments = roots[rows[pre_indexes]]
        weights = vals[pre_indexes]
        reversal = reversal_map[syn_model]
        post_synapses = get_synapses(post_compartment)
        
        for (i, (pre_comp, weight)) in enumerate(zip(pre_compartments, weights))
            model_copy = remake(syn_model, subscriptions = [pre_comp], gbar = weight,
                                name = Symbol(nameof(syn_model), "_$i"))
            push!(post_synapses, Synapse(model_copy, reversal))
            push!(model_copies, model_copy)
        end

        post_compartment = remake(post_compartment; synapses = post_synapses)

        # Directly connect Vₘ between neurons; must happen post hoc

        for (model_copy, pre_comp) in zip(model_copies, pre_compartments)
            vars = namespace_variables(getproperty(post_compartment, nameof(model_copy)))
            Vx = find_voltage(vars, isextrinsic)
            push!(voltage_fwds, pre.Vₘ ~ Vx)
            push!(defs, Vx => pre.Vₘ)
        end

        new_compartments[i] = post_compartment
    end
    union!(eqs, voltage_fwds)
    return new_compartments
end

"""
$(TYPEDSIGNATURES)

Basic constructor for a `NeuronalNetworkSystem`.

"""
function NeuronalNetworkSystem(topology::NetworkTopology, reversal_map,
                               extensions::Vector{<:AbstractTimeDependentSystem} = AbstractTimeDependentSystem[];
                               defaults = Dict(),
                               name::Symbol = Base.gensym(:Network)) where {T}

    reversal_map = reversal_map isa AbstractDict ? reversal_map : Dict(reversal_map)

    gen = GeneratedCollections()
    (; eqs, dvs, ps, systems, observed, defs) = gen
    comps = compartments(topology) 
    multigraph = graph(topology)
    for synaptic_model in synaptic_models(topology)
        comps = connect_synapses!(gen, synaptic_model, comps, multigraph)
    end

    union!(systems, extensions, compartments)
    return NeuronalNetworkSystem(eqs, t, collect(dvs), collect(ps), observed, name, systems,
                                 defaults, topology, reversal_map, extensions;
                                 checks = false)
end

get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)
get_topology(x::AbstractNeuronalNetworkSystem) = getfield(x, :topology)

Base.eltype(ns::NeuronalNetworkSystem) = AbstractCompartmentSystem
Base.length(ns::NeuronalNetworkSystem) = length(neurons(get_topology(ns)))

function Base.iterate(ns::NeuronalNetworkSystem, state = 1)
    state > length(ns) && return nothing
    return (neurons(get_topology(ns))[state], state + 1)
end

function Base.iterate(rev_ns::Iterators.Reverse{NeuronalNetworkSystem},
                      state = length(rev_ns.itr))
    state < 1 && return nothing
    ns = rev_ns.itr
    return (neurons(get_topology(ns))[state], state - 1)
end

function Base.getindex(ns::NeuronalNetworkSystem, i)
    return neurons(get_topology(ns))[i]
end

Base.firstindex(ns::NeuronalNetworkSystem) = 1
Base.lastindex(ns::NeuronalNetworkSystem) = length(get_topology(ns))

function Base.convert(::Type{ODESystem}, nnsys::NeuronalNetworkSystem)
    dvs = get_states(nnsys)
    ps = get_ps(nnsys)
    eqs = get_eqs(nnsys)
    defs = get_defaults(nnsys)
    ccbs = get_continuous_events(nnsys)
    dcbs = get_discrete_events(nnsys)
    syss = convert.(ODESystem, get_systems(nnsys))
    return ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(nnsys), systems = syss,
                     discrete_events = dcbs, continuous_events = ccbs,
                     checks = CheckComponents)
end
