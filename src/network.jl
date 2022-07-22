
struct NetworkTopology
    multigraph::Dict{ConductanceSystem, SparseMatrixCSC{Num,Int64}}
    neuron_map::Dict{AbstractCompartmentSystem, UnitRange{Int64}}
    compartments::Vector{CompartmentSystem}
    denamespaced::Vector{Symbol}
end

function NetworkTopology(neurons::Vector{<:AbstractCompartmentSystem},
                         synaptic_models::Vector{<:AbstractConductanceSystem})
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
    multigraph = Dict([x => sparse(Int64[], Int64[], Num[], n, n) for x in synaptic_models])

    return NetworkTopology(multigraph, neuron_map, all_compartments, denamespaced)
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
    return NetworkTopology(multigraph, neuron_map, all_compartments, denamespaced)
end

neurons(topology::NetworkTopology) = [keys(topology.neuron_map)...]
vertices(topology::NetworkTopology) = getfield(topology, :compartments)
graph(topology::NetworkTopology) = getfield(topology, :multigraph)

function add_synapse!(topology, pre, post, synaptic_model)
    src = find_compsys(pre, topology)
    dst = find_compsys(post, topology)
    g = graph(topology)[synaptic_model]
    g[src, dst] = get_gbar(synaptic_model)
end

function remove_synapse!(topology, pre, post, synaptic_model)
    src = find_compsys(pre, topology)
    dst = find_compsys(post, topology)
    g = graph(topology)[synaptic_model]
    g[src,dst] = zero(Num)
    dropzeros!(g)
end

function add_layer!(topology, synaptic_model, g = SimpleDiGraph(length(vertices(topology))), default_weight = get_gbar(synaptic_model))
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
    states::Vector{Num}
    ps::Vector{Num}
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
                                   topology, reversal_map, extensions; checks = false)
        if checks
            # placeholder
        end
        new(eqs, iv, states, ps, observed, name, systems, defaults, topology, reversal_map, extensions)
    end
end

"""
$(TYPEDSIGNATURES)

Basic constructor for a `NeuronalNetworkSystem`.

"""
function NeuronalNetworkSystem(
    topology, reversal_map,
    extensions::Vector{<:AbstractTimeDependentSystem} = AbstractTimeDependentSystem[];
    defaults = Dict(), name::Symbol = Base.gensym(:Network))
    
    reversal_map = reversal_map isa AbstractDict ? reversal_map : Dict(reversal_map)
    eqs, dvs, ps, observed = Equation[], Num[], Num[], Equation[]
    voltage_fwds = Set{Equation}()
    multigraph = graph(topology)
    compartments = vertices(topology)

    for synaptic_class in keys(multigraph)
        g = multigraph[synaptic_class]
        rows = rowvals(g) # row numbers for each stored value
        vals = nonzeros(g) # stored values
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
                Vxs = filter(x -> isvoltage(x) && isextrinsic(x),
                             MTK.namespace_variables(getproperty(post_compartment, nameof(new_synaptic_class))))
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
                #####
                post_dynamics = get_dynamics(post_compartment)
                new_dynamics = @set post_dynamics.synaptic_channels = class_copies
                @set! new_dynamics.synaptic_reversals = new_revs
                post_compartment = SciMLBase.remake(post_compartment, dynamics = new_dynamics)
                for (class_copy, pre) in zip(class_copies, pre_compartments)
                    Vx = getproperty(post_compartment, nameof(class_copy)).Vₓ
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
        idx_range = topology.neuron_map[neuron]
        if neuron isa MultiCompartmentSystem
            mctop = get_topology(neuron)
            # subcompartments can't be namespaced!
            @set! mctop.compartments = MTK.rename.(compartments[idx_range],
                                                   topology.denamespaced[idx_range])
            neuron = MultiCompartmentSystem(neuron, topology = mctop)
        else
            neuron = only(compartments[idx_range])
        end
        newmap[neuron] = idx_range
    end
    systems::Vector{AbstractTimeDependentSystem} = union(extensions, collect(keys(newmap)))
    return NeuronalNetworkSystem(eqs, t, dvs, ps, observed, name, systems, defaults,
                                 topology, reversal_map, extensions; checks = false)
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
#
#function Base.setindex!(ns::NeuronalNetworkSystem, cond::ConductanceSystem, pre::T1, post::T2) where {T1<:AbstractCompartmentSystem, T2<:AbstractCompartmentSystem}
#    add_synapse!(ns, pre, post, cond)
#end

function Base.convert(::Type{ODESystem}, nnsys::NeuronalNetworkSystem)
    dvs = get_states(nnsys)
    ps  = get_ps(nnsys)
    eqs = get_eqs(nnsys)
    defs = get_defaults(nnsys)
    syss = convert.(ODESystem, get_systems(nnsys))
    return ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(nnsys), systems = syss)
end


