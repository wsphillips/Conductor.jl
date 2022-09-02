
struct NeuronalMultigraph
    multigraph::Dict{Any, SparseMatrixCSC{Num,Int64}}
    neuron_map::Dict{AbstractCompartmentSystem, UnitRange{Int64}}
    denamespaced::Vector{Symbol}
end

struct NetworkTopology{T,F}
    graph::T
    compartments::Vector{CompartmentSystem{F}}
end

function NetworkTopology(neurons::Vector{CompartmentSystem{LIF}}, graph = SimpleDiGraph(lenght(neurons)); default_weight = 1.0)
    n = length(neurons)
    graph = adjacency_matrix(graph)*default_weight
    NetworkTopology(graph, neurons)
end

function NetworkTopology(neurons::Vector{<:AbstractCompartmentSystem},
                         synaptic_models::Vector{<:AbstractConductanceSystem})

    neuron_map = Dict{AbstractCompartmentSystem, UnitRange{Int64}}()
    start = 1
    all_compartments = []
    denamespaced = Symbol[]

    for neuron in neurons
        comps = compartments(neuron)
        append!(all_compartments, comps)
        denamespaced_comps = nameof.(compartments(neuron, namespace = false))
        append!(denamespaced, denamespaced_comps)
        n_comps = length(comps)
        neuron_map[neuron] = start:(start + n_comps - 1)
        start += n_comps
    end

    n = length(all_compartments)
    multigraph = Dict([x => sparse(Int64[], Int64[], Num[], n, n) for x in synaptic_models])
    graph = NeuronalMultigraph(multigraph, neuron_map, denamespaced)
    return NetworkTopology(graph, [all_compartments...])
end

function NetworkTopology(g::SimpleDiGraph, neurons::Vector{<:AbstractCompartmentSystem},
        synaptic_model::AbstractConductanceSystem, default_weight = get_gbar(synaptic_model))
    
    neuron_map = Dict{AbstractCompartmentSystem, UnitRange{Int64}}()
    start = 1
    all_compartments = CompartmentSystem[]
    denamespaced = Symbol[]

    for neuron in neurons
        comps = compartments(neuron)
        append!(all_compartments, comps)
        denamespaced_comps = nameof.(compartments(neuron, namespace = false))
        append!(denamespaced, denamespaced_comps)
        n_comps = length(comps)
        neuron_map[neuron] = start:(start + n_comps - 1)
        start += n_comps
    end

    n = length(all_compartments)
    @assert length(vertices(g)) == length(neurons)
    multigraph = Dict(synaptic_model => adjacency_matrix(g)*default_weight)
    graph = NeuronalMultigraph(multigraph, neuron_map, denamespaced)

    return NetworkTopology(graph, all_compartments)
end

function neurons(topology::NetworkTopology{NeuronalMultigraph,<:CompartmentForm})
    [keys(neuron_map(topology))...]
end

denamespaced(topology::NetworkTopology{NeuronalMultigraph}) = getfield(topology,:graph).denamespaced
neuron_map(topology::NetworkTopology{NeuronalMultigraph}) = getfield(topology,:graph).neuron_map
neurons(topology::NetworkTopology{T,F}) where {T,F} = vertices(topology)
vertices(topology::NetworkTopology) = getfield(topology, :compartments)
graph(topology::NetworkTopology) = getfield(topology, :graph)
graph(topology::NetworkTopology{NeuronalMultigraph}) = getfield(topology, :graph).multigraph
function add_synapse!(topology, pre, post, synaptic_model::ConductanceSystem)
    src = find_compsys(pre, topology)
    dst = find_compsys(post, topology)
    g = graph(topology)[synaptic_model]
    g[src, dst] = get_gbar(synaptic_model)
end

function add_synapse!(topology, pre::AbstractCompartmentSystem, post::AbstractCompartmentSystem, weight::Real)
    src = find_compsys(pre, topology)
    dst = find_compsys(post, topology)
    graph(topology)[src, dst] = weight
end

function add_synapse!(topology, pre::Int64, post::Int64, weight::Real)
    graph(topology)[pre, post] = weight
end

function remove_synapse!(topology, pre, post, synaptic_model)
    src = find_compsys(pre, topology)
    dst = find_compsys(post, topology)
    g = graph(topology)[synaptic_model]
    g[src,dst] = zero(Num)
    dropzeros!(g)
end

function add_layer!(topology, synaptic_model::ConductanceSystem, g = SimpleDiGraph(length(vertices(topology))), default_weight = get_gbar(synaptic_model))
    push!(graph(topology), synaptic_model => adjacency_matrix(g)*default_weight)
end

Base.eltype(nt::NetworkTopology) = AbstractCompartmentSystem
Base.length(nt::NetworkTopology) = length(neurons(nt))

function Base.iterate(nt::NetworkTopology, state=1)
    state > length(nt) && return nothing
    return (neurons(nt)[state], state+1)
end

function Base.iterate(rev_nt::Iterators.Reverse{NetworkTopology}, state=length(rev_nt.itr))
    state < 1 && return nothing
    nt = rev_nt.itr
    return (neurons(nt)[state], state-1)
end

function Base.getindex(nt::NetworkTopology, i)
    return neurons(nt)[i]
end

Base.firstindex(nt::NetworkTopology) = 1
Base.lastindex(nt::NetworkTopology) = length(nt)

function Base.setindex!(nt::NetworkTopology, g::SimpleDiGraph, cond::ConductanceSystem)
    add_layer!(nt, cond, g)
end

function Base.setindex!(nt::NetworkTopology, cond::ConductanceSystem, pre::T1, post::T2) where {T1<:AbstractCompartmentSystem, T2<:AbstractCompartmentSystem}
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
    continuous_events::Vector{MTK.SymbolicContinuousCallback}
    discrete_events::Vector{MTK.SymbolicDiscreteCallback}
    # Conductor fields
    topology::NetworkTopology
    reversal_map::Dict
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    function NeuronalNetworkSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                                   continuous_events, discrete_events, topology, reversal_map,
                                   extensions; checks = false)
        if checks
            # placeholder
        end
        new(eqs, iv, states, ps, observed, name, systems, defaults, continuous_events,
            discrete_events, topology, reversal_map, extensions)
    end
end

"""
$(TYPEDSIGNATURES)

Basic constructor for a `NeuronalNetworkSystem`.

"""
function NeuronalNetworkSystem(topology::NetworkTopology{T, HodgkinHuxley}, reversal_map,
    extensions::Vector{<:AbstractTimeDependentSystem} = AbstractTimeDependentSystem[];
    defaults = Dict(), name::Symbol = Base.gensym(:Network)) where T
    
    reversal_map = reversal_map isa AbstractDict ? reversal_map : Dict(reversal_map)
    eqs, dvs, ps, observed = Equation[], Num[], Num[], Equation[]
    ccbs, dcbs = MTK.SymbolicContinuousCallback[], MTK.SymbolicDiscreteCallback[]
    voltage_fwds = Set{Equation}()
    multigraph = graph(topology)
    compartments = vertices(topology)

    for synaptic_class in keys(multigraph)
        g = multigraph[synaptic_class]
        rows = rowvals(g)
        vals = nonzeros(g)
        for i in axes(g, 2) # for each set of presynaptic neurons per postsynaptic neuron
            pre_indexes = nzrange(g, i)
            iszero(length(pre_indexes)) && continue
            post_compartment = compartments[i]
            pre_compartments = compartments[rows[pre_indexes]] # use view?
            weights = vals[pre_indexes] # use view?
            new_revs = union(get_synaptic_reversals(post_compartment),
                             reversal_map[synaptic_class])
            if isaggregate(synaptic_class)
                new_synaptic_class = ConductanceSystem(synaptic_class,
                                                       subscriptions = pre_compartments) 
                post_dynamics = get_dynamics(post_compartment)
                new_dynamics = @set post_dynamics.synaptic_channels = [new_synaptic_class]
                @set! new_dynamics.synaptic_reversals = new_revs
                post_compartment = SciMLBase.remake(post_compartment, dynamics = new_dynamics)
                vars = MTK.namespace_variables(getproperty(post_compartment,
                                                           nameof(new_synaptic_class)))
                Vxs = find_voltage(vars, isextrinsic)
                for (Vx, pre) in zip(Vxs, pre_compartments)
                    push!(voltage_fwds, pre.Vₘ ~ Vx)
                    push!(defaults, Vx => pre.Vₘ)
                end
            else
                class_copies = [] # clone the synapse model for each presynaptic compartment
                for (x,y) in zip(pre_compartments, weights)
                    class_copy = ConductanceSystem(synaptic_class, subscriptions = [x],
                                                   gbar = y,
                                                   name = namegen(nameof(synaptic_class)),
                                                   defaults = Dict(y => getdefault(y)))
                    push!(class_copies, class_copy)
                end
                post_dynamics = get_dynamics(post_compartment)
                new_dynamics = @set post_dynamics.synaptic_channels = class_copies
                @set! new_dynamics.synaptic_reversals = new_revs
                post_compartment = SciMLBase.remake(post_compartment, dynamics = new_dynamics)
                for (class_copy, pre) in zip(class_copies, pre_compartments)
                    vars = MTK.namespace_variables(getproperty(post_compartment, nameof(class_copy)))
                    Vx = find_voltage(vars, isextrinsic) 
                    push!(voltage_fwds, pre.Vₘ ~ Vx)
                    push!(defaults, Vx => pre.Vₘ)
                end
            end
            compartments[i] = post_compartment
        end
    end

    union!(eqs, voltage_fwds)

    newmap = Dict{AbstractCompartmentSystem, UnitRange{Int64}}()
    # Construct revised neurons
    for neuron in neurons(topology)
        idx_range = neuron_map(topology)[neuron]
        if neuron isa MultiCompartmentSystem
            mctop = get_topology(neuron)
            # subcompartments can't be namespaced!
            @set! mctop.compartments = MTK.rename.(compartments[idx_range],
                                                   denamespaced(topology)[idx_range])
            neuron = MultiCompartmentSystem(neuron, topology = mctop)
        else
            neuron = only(compartments[idx_range])
        end
        newmap[neuron] = idx_range
    end
    systems::Vector{AbstractTimeDependentSystem} = union(extensions, collect(keys(newmap)))
    return NeuronalNetworkSystem(eqs, t, dvs, ps, observed, name, systems, defaults,
                                 ccbs, dcbs, topology, reversal_map, extensions; checks = false)
end

function NeuronalNetworkSystem(topology::NetworkTopology{T, LIF}, 
        extensions::Vector{<:AbstractTimeDependentSystem} = AbstractTimeDependentSystem[];
        defaults = Dict(), name::Symbol = Base.gensym(:Network)) where T

    g = graph(topology) 
    neurons = vertices(topology)
    n = length(neurons)
    @parameters W[1:n,1:n] = g 
    #@variables (Syn_inp(t))[1:n] = zeros(n)
    reversal_map = Dict() 
    gen = GeneratedCollections(systems = union(neurons, extensions),
                               ps = Set(scalarize(W)))

    (; eqs, dvs, ps, systems, observed) = gen
    
    affect_eqs = Equation[]
    spike_checks = [neuron.V >= neuron.V_th for neuron in neurons]

    for (i, neuron) in enumerate(neurons)
        push!(affect_eqs, neuron.I ~ neuron.I + sum(scalarize(W[:,i] .* spike_checks)))
    end
    # equations _are_ order-dependent here. We can do this more nicely with a generic func affect 
    for (i, neuron) in enumerate(neurons)
        push!(affect_eqs, neuron.V ~ ifelse(spike_checks[i], neuron.V_rest, neuron.V))
    end

    condition = reduce(|, spike_checks) # also tried with periodic callbacks--same issue
    dcbs = [MTK.SymbolicDiscreteCallback(condition, affect_eqs)]
    ccbs = MTK.SymbolicContinuousCallback[]
    return NeuronalNetworkSystem(eqs, t, collect(dvs), collect(ps), observed, name, systems,
                                 defaults, ccbs, dcbs, topology, reversal_map, extensions;
                                 checks = false)
end

get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)
get_topology(x::AbstractNeuronalNetworkSystem) = getfield(x, :topology)

Base.eltype(ns::NeuronalNetworkSystem) = AbstractCompartmentSystem
Base.length(ns::NeuronalNetworkSystem) = length(neurons(get_topology(ns)))

function Base.iterate(ns::NeuronalNetworkSystem, state=1)
    state > length(ns) && return nothing
    return (neurons(get_topology(ns))[state], state+1)
end

function Base.iterate(rev_ns::Iterators.Reverse{NeuronalNetworkSystem}, state=length(rev_ns.itr))
    state < 1 && return nothing
    ns = rev_ns.itr
    return (neurons(get_topology(ns))[state], state-1)
end

function Base.getindex(ns::NeuronalNetworkSystem, i)
    return neurons(get_topology(ns))[i]
end

Base.firstindex(ns::NeuronalNetworkSystem) = 1
Base.lastindex(ns::NeuronalNetworkSystem) = length(get_topology(ns))

#function Base.setindex!(ns::NeuronalNetworkSystem, g::SimpleDiGraph, cond::ConductanceSystem)
#    add_layer!(ns, cond, g)
#end
#function Base.setindex!(ns::NeuronalNetworkSystem, cond::ConductanceSystem, pre::T1, post::T2) where {T1<:AbstractCompartmentSystem, T2<:AbstractCompartmentSystem}
#    add_synapse!(ns, pre, post, cond)
#end

function Base.convert(::Type{ODESystem}, nnsys::NeuronalNetworkSystem)
    dvs = get_states(nnsys)
    ps  = get_ps(nnsys)
    eqs = get_eqs(nnsys)
    defs = get_defaults(nnsys)
    ccbs = get_continuous_events(nnsys)
    dcbs = get_discrete_events(nnsys)
    syss = convert.(ODESystem, get_systems(nnsys))
    return ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(nnsys), systems = syss, discrete_events = dcbs, continuous_events = ccbs, checks = CheckComponents)
end

