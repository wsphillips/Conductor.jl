
struct NetworkTopology
    multigraph::Dict{ConductanceSystem, SparseMatrixCSC{Num,Int64}}
    neuron_map::Dict{Any, UnitRange{Int64}}
    compartments::Vector{CompartmentSystem}
end

function NetworkTopology(neurons::Vector{<:AbstractCompartmentSystem},
                         synaptic_models::Vector{<:AbstractConductanceSystem})
    neuron_map = Dict()
    start = 1
    all_compartments = []
    for neuron in neurons
        comps = compartments(neuron)
        append!(all_compartments, comps)
        n_comps = length(comps)
        neuron_map[neuron] = start:(start + n_comps - 1)
        start += n_comps
    end
    n = length(all_compartments)
    multigraph = Dict([x => sparse(Int64[], Int64[], Num[], n, n) for x in synaptic_models])
    return NetworkTopology(multigraph, neuron_map, all_compartments)
end

function NetworkTopology(g::SimpleDiGraph, neurons::Vector{<:AbstractCompartmentSystem},
        synaptic_model::AbstractConductanceSystem, default_weight = get_gbar(synaptic_model))
    compartments = CompartmentSystem[]

    for neuron in neurons
         # FIXME!!! need list of compartments for graph and index mapping to reconstruct neurons 
    end
    @assert length(vertices(g)) == length(neurons)
    multigraph = Dict(synaptic_model => adjacency_matrix(g)*default_weight)
    return NetworkTopology(multigraph, neurons)
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

    eqs = Equation[]
    dvs = Set{Num}()
    ps  = Set{Num}()
    observed = Equation[]
    voltage_fwds = Set{Equation}()
    multigraph = graph(topology)
    compartments = vertices(topology)

    for synaptic_class in keys(multigraph)
        g = multigraph[synaptic_class]
        rows = rowvals(g) # row numbers for each stored value
        vals = nonzeros(g) # stored values
        for i in axes(g, 2) # for each set of presynaptic neurons per postsynaptic neuron
            post_compartment = compartments[i]
            pre_compartments = compartments[rows[nzrange(g, i)]]
            pre_compartments = pre_compartments isa Vector ? pre_compartments :
                                                             [pre_compartments]
            weights = vals[nzrange(g,i)]
            new_revs = union(get_synaptic_reversals(post_compartment),
                             reversal_map[synaptic_class])
            post_compartment = CompartmentSystem(post_compartment,
                                                 synaptic_reversals = new_revs)
            if isaggregate(synaptic_class)
                synaptic_class = ConductanceSystem(synaptic_class,
                                                   subscriptions = pre_compartments) 
                post_compartment = CompartmentSystem(post_compartment,
                                                     synaptic_channels = [synaptic_class])
                Vxs = filter(x -> isvoltage(x) && isextrinsic(x),
                             namespace_variables(post_compartment.synaptic_class))
                for (Vx, pre) in zip(Vxs, pre_compartments)
                    push!(voltage_fwds, pre.Vₘ ~ Vx)
                    push!(defaults, Vx => pre.Vₘ)
                end
                compartments[i] = post_compartment
            else
                class_copies = [] # clone the synapse model for each presynaptic compartment
                for (x,y) in zip(pre_compartments, weights)
                    class_copy = ConductanceSystem(synaptic_class, subscriptions = [x],
                                                   name = namegen(nameof(synaptic_class)),
                                                   defaults = Dict(y => getdefault(y)))
                    push!(class_copies, class_copy)
                end
                post_compartment = CompartmentSystem(post_compartment,
                                                     synaptic_channels = class_copies)
                for (class_copy, pre) in zip(class_copies, pre_compartments)
                    Vx = getproperty(post_compartment, nameof(class_copy)).Vₓ
                    push!(voltage_fwds, pre.Vₘ ~ Vx)
                    push!(defaults, Vx => pre.Vₘ)
                end
                compartments[i] = post_compartment
            end
        end
    end
    union!(eqs, voltage_fwds)

    newmap = Dict()
    # reconstruct the multicompartment neurons
    for neuron in neurons(topology)
        idx_range = topology.neuron_map[neuron]
        if neuron isa MultiCompartmentSystem
            mctop = get_topology(neuron)
            @set! mctop.compartments = compartments[idx_range]
            neuron = MultiCompartmentSystem(neuron, topology = mctop)
        else
            neuron = only(compartments[idx_range])
        end
        newmap[neuron] = idx_range
    end
    topology = NetworkTopology(topology.multigraph, newmap, compartments)
    systems::Vector{AbstractTimeDependentSystem} = union(extensions, neurons(topology))

    return NeuronalNetworkSystem(eqs, t, collect(dvs), collect(ps), observed, name, systems, defaults,
                                   topology, reversal_map, extensions; checks = false)
end

get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)
get_topology(x::AbstractNeuronalNetworkSystem) = getfield(x, :topology)

function Base.convert(::Type{ODESystem}, nnsys::NeuronalNetworkSystem)
    dvs = get_states(nnsys)
    ps  = get_ps(nnsys)
    eqs = get_eqs(nnsys)
    defs = get_defaults(nnsys)
    syss = convert.(ODESystem, collect(get_systems(nnsys)))
    return ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(nnsys), systems = syss)
end

#=
function Base.:(==)(sys1::NeuronalNetworkSystem, sys2::NeuronalNetworkSystem)
    sys1 === sys2 && return true
    iv1 = get_iv(sys1)
    iv2 = get_iv(sys2)
    isequal(iv1, iv2) &&
    isequal(nameof(sys1), nameof(sys2)) &&
    _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
    _eq_unordered(get_states(sys1), get_states(sys2)) &&
    _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
    all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2)))
end
=#

function reset_synapses!(x::CompartmentSystem)
    empty!(get_synapses(x))
    empty!(get_synaptic_reversals(x))
    return nothing
end

function reset_synapses!(x::MultiCompartmentSystem)
    foreach(reset_synapses!, get_compartments(x))
end


