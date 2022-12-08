
struct NetworkTopology
    multigraph::Dict{Any, SparseMatrixCSC{Float64, Int64}}
    neurons::Vector
end

synaptic_systems(topology::NetworkTopology) = keys(graph(topology))

function NetworkTopology(neurons::Vector, synaptic_systems::Vector)
    m = length(neurons)
    n = sum(length(neuron) for neuron in neurons)
    multigraph = Dict(x => sparse(Int64[], Int64[], Float64[], n, m) for x in synaptic_systems)
    return NetworkTopology(multigraph, neurons)
end

function NetworkTopology(g::SimpleDiGraph, neurons, synaptic_system,
                         default_weight = get_gbar(synaptic_system))
    multigraph = Dict(synaptic_system => adjacency_matrix(g) * default_weight)
    return NetworkTopology(multigraph, neurons)
end

neurons(topology::NetworkTopology) = getfield(topology, :neurons)

# returns the _individual_ compartments as a single concatenated vector
function compartments(topology::NetworkTopology; namespace = true)
    foldl(vcat, compartments(neuron; namespace) for neuron in neurons(topology))
end

root_compartments(topology::NetworkTopology) = [first(neuron) for neuron in neurons(topology)]
graph(topology::NetworkTopology) = getfield(topology, :multigraph)

find_source(pre, topology) = findfirst(isequal(first(pre)), root_compartments(topology))
find_target(post, topology) = findfirst(isequal(post), compartments(topology))

function add_synapse!(topology, pre, post, synaptic_system, weight)
    src = find_source(pre, topology)
    dst = find_target(post, topology)
    g = graph(topology)[synaptic_system]
    g[dst, src] = weight # underlying matrix is transposed
end

# FIXME: need updates
#=
function remove_synapse!(topology, pre, post, synaptic_system)
    src = find_source(pre, topology)
    dst = find_target(post, topology)
    g = graph(topology)[synaptic_system]
    g[dst, src] = zero(Num)
    dropzeros!(g)
end

function add_layer!(topology, synaptic_system::ConductanceSystem,
                    g = SimpleDiGraph(length(vertices(topology))),
                    default_weight = get_gbar(synaptic_system))
    push!(graph(topology), synaptic_system => adjacency_matrix(g) * default_weight)
end
=#

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

function Base.getindex(nt::NetworkTopology, i::Int64)
    return neurons(nt)[i]
end

function Base.getindex(nt::NetworkTopology, sys::AbstractSystem)
    return nt.multigraph[sys]
end

Base.firstindex(nt::NetworkTopology) = 1
Base.lastindex(nt::NetworkTopology) = length(nt)

function Base.setindex!(nt::NetworkTopology, g::SimpleDiGraph, cond::ConductanceSystem)
    add_layer!(nt, cond, g)
end

function Base.setindex!(nt::NetworkTopology, cond::ConductanceSystem, pre::T1,
                        post::T2) where {T1 <: AbstractCompartmentSystem,
                                         T2 <: AbstractCompartmentSystem}
    gbar = get_gbar(cond)
    add_synapse!(nt, pre, post, cond, getdefault(gbar))
end

function Base.setindex!(nt::NetworkTopology, cw::Tuple{ConductanceSystem, Any}, pre, post)
    cond, new_weight = cw
    maxval = getdefault(get_gbar(cond))
    val = new_weight/maxval
    add_synapse!(nt, pre, post, cond, val)
end

function Base.setindex!(nt::NetworkTopology, cond::ConductanceSystem, pre, post)
    add_synapse!(nt, pre, post, cond, 1.0)
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

synaptic_systems(sys::NeuronalNetworkSystem) = synaptic_systems(get_topology(sys))
compartments(sys::NeuronalNetworkSystem) = compartments(get_topology(sys))

function connect_synapses!(gen, syn_model, comps, topology, reversal_map)
        # this method assumes weights scale the size of the event/alpha)
    new_compartments = deepcopy(comps)
    reversal = reversal_map[syn_model]
    for (i,comp) in enumerate(new_compartments)
        post_synapses = get_synapses(comp)
        push!(post_synapses, Synapse(syn_model, reversal))
        new_compartments[i] = remake(comp; synapses = post_synapses)
    end
    return new_compartments
end

function connect_synapses!(gen, syn_model::Union{CurrentSystem{T},ConductanceSystem{T}},
                           comps, topology, reversal_map) where {T<:IntegratedSynapse}

    (; eqs, systems, defs) = gen
    new_compartments = deepcopy(comps)
    roots = root_compartments(topology)
    g = permutedims(graph(topology)[syn_model]) # easier to do with rows = pre, cols = post
    rows = rowvals(g)
    vals = nonzeros(g)
    voltage_fwds = Set{Equation}()
    @parameters W 
    # Synaptic conductance models are duplicated here and individually connected, relying
    # on users to specify extrinsic voltage symbols. In the future, when symbolic arrays and
    # component types become available we will probably rewrite how this works.
    for i in axes(g, 2) # for each set of presynaptic neurons per postsynaptic neuron
        pre_indexes = rows[nzrange(g, i)]
        isempty(pre_indexes) && continue # skip neurons with no synaptic input

        model_copies = []
        post_compartment = new_compartments[i]
        pre_compartments = roots[rows[pre_indexes]]
        weights = vals[pre_indexes]
        reversal = reversal_map[syn_model]
        post_synapses = get_synapses(post_compartment) # currently, this gets mutated
        
        # duplicate to allow for multiple presynaptic inputs of the same type
        for (i, (pre_comp, weight)) in enumerate(zip(pre_compartments, weights))
            model_copy = remake(syn_model, name = Symbol(nameof(syn_model), "_$i"),
                                defaults = Dict(W => weight))
            push!(post_synapses, Synapse(model_copy, reversal))
            push!(model_copies, model_copy)
        end

        post_compartment = remake(post_compartment; synapses = post_synapses)

        # Directly connect Vₘ between neurons; must happen post hoc

        for (mcopy, pre_comp) in zip(model_copies, pre_compartments)
            vars = namespace_variables(getproperty(post_compartment, nameof(mcopy)))
            Vx = find_voltage(vars, isextrinsic)
            push!(voltage_fwds, pre_comp.Vₘ ~ Vx)
            push!(defs, Vx => pre_comp.Vₘ)
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
                               name::Symbol = Base.gensym(:Network))

    reversal_map = reversal_map isa AbstractDict ? reversal_map : Dict(reversal_map)

    gen = GeneratedCollections()
    (; eqs, dvs, ps, systems, observed, defs) = gen
    comps = compartments(topology; namespace = false)
    for sys in synaptic_systems(topology)
        comps = connect_synapses!(gen, sys, comps, topology, reversal_map)
    end

    for neuron in neurons(topology)
        if typeof(neuron) <: MultiCompartmentSystem
            mctop = get_topology(neuron)
            @set! mctop.compartments = comps[1:length(neuron)]
            push!(systems, remake(neuron, topology = mctop))
            length(comps) == length(neuron) && break
            comps = comps[length(neuron)+1:end]
        else
            push!(systems, popfirst!(comps))
        end
    end
    union!(systems, extensions)
    merge!(defs, defaults)
    return NeuronalNetworkSystem(eqs, t, collect(dvs), collect(ps), observed, name, systems,
                                 defs, topology, reversal_map, extensions;
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
