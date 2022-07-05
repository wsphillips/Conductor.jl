
abstract type AbstractSynapse end # <: AbstractEdge{T} end

"""
$(TYPEDEF)

A synapse (edge) between two neurons in a network.

$(FIELDS)
"""
struct Synapse <: AbstractSynapse
    "Presynaptic neuron/subcompartment."
    source::CompartmentSystem
    "Postsynaptic neuron/subcompartment."
    target::CompartmentSystem
    "Synaptic conductance model of the synapse."
    conductance::AbstractConductanceSystem
    "Synaptic equilibrium/reversal potential."
    reversal::Num
end

"""
$(TYPEDSIGNATURES)

Basic constructor for a directed `Synapse` (edge) between two compartments.

Specified as a pair of compartments:

`presynaptic::CompartmentSystem => postsynaptic::CompartmentSystem`.
"""
function Synapse(pre_to_post::Pair, conductance, reversal)
    return Synapse(pre_to_post.first, pre_to_post.second, conductance, reversal) 
end

presynaptic(x::Synapse) = getfield(x, :source)
postsynaptic(x::Synapse) = getfield(x, :target)
class(x::Synapse) = getfield(x, :conductance)
reversal(x::Synapse) = getfield(x, :reversal)

struct NetworkTopology
    multigraph::Dict{ConductanceSystem, SparseMatrixCSC{Num,Int64}}
    neurons::Vector{AbstractCompartmentSystem}
end

vertices(topology::NetworkTopology) = getfield(topology, :neurons)
graph(topology::NetworkTopology) = getfield(topology, :multigraph)

# create an topology with no connections
function NetworkTopology(neurons::Vector{<:AbstractCompartmentSystem},
                         synaptic_models::Vector{<:AbstractConductanceSystem})
    n = length(neurons)
    multigraph = Dict([x => sparse(Int64[], Int64[], Num[], n, n) for x in synaptic_models])
    return NetworkTopology(multigraph, neurons)
end

# create a pre-specified topology (from a graph template)
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
    extensions::Vector{AbstractTimeDependentSystem}
    function NeuronalNetworkSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                                   topology, reversal_map; checks = false)
        if checks
            # placeholder
        end
        new(eqs, iv, states, ps, observed, name, systems, defaults, topology, reversal_map)
    end
end

"""
$(TYPEDSIGNATURES)

Basic constructor for a `NeuronalNetworkSystem`.

"""
function NeuronalNetworkSystem(
    topology, reversal_map,
    extensions::Vector{<:AbstractTimeDependentSystem} = AbstractTimeDependentSystem[];
    defs = Dict(), name::Symbol = Base.gensym(:Network))

    eqs = Equation[]
    dvs = Set{Num}()
    ps  = Set{}()
    observed = Equation[]
    systems = Vector{AbstractTimeDependentSystem}()
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
                                                     synapses = [synaptic_class])
                Vxs = filter(x -> isvoltage(x) && isextrinsic(x),
                             namespace_variables(post_compartment.synaptic_class))
                for (Vx, pre) in zip(Vxs, pre_compartments)
                    push!(voltage_fwds, pre.Vₘ ~ Vx)
                    push!(defs, Vx => pre.Vₘ)
                end
            else
                class_copies = [] # clone the synapse model for each presynaptic compartment
                for (x,y) in zip(pre_compartments, weights)
                    class_copy = ConductanceSystem(synaptic_class, subscriptions = [x],
                                                   name = namegen(nameof(synaptic_class)),
                                                   defaults = Dict(y => getdefault(y)))
                    push!(class_copies, class_copy)
                end
                post_compartment = CompartmentSystem(post_compartment,
                                                     synapses = class_copies)
                for (class_copy, pre) in zip(class_copies, pre_compartments)
                    Vx = post_compartment.class_copy.Vₓ
                    push!(voltage_fwds, pre.Vₘ ~ Vx)
                    push!(defs, Vx => pre.Vₘ)
                end
            end
            compartments[i] = post_compartment
        end
    end
    union!(eqs, voltage_fwds)
    # FIXME: must go back and re-run MultiCompartment constructors since component
    # compartments will be new/updated
    return NeuronalNetworkSystem(eqs, t, states, ps, observed, name, systems, defaults,
                                   topology, reversal_map; checks = false)
end

get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)
get_topology(x::AbstractNeuronalNetworkSystem) = getfield(x, :topology)

function Base.convert(::Type{ODESystem}, nnsys::NeuronalNetworkSystem)
    states, params, eqs, defs, allneurons = build_toplevel(nnsys)
    all_systems = map(x -> convert(ODESystem, x), allneurons)
    odesys = ODESystem(eqs, t, states, params; defaults = defs, name = nameof(nnsys))
    return compose(odesys, all_systems)
end

function reset_synapses!(x::CompartmentSystem)
    empty!(get_synapses(x))
    empty!(get_synaptic_reversals(x))
    return nothing
end

function reset_synapses!(x::MultiCompartmentSystem)
    foreach(reset_synapses!, get_compartments(x))
end

function build_toplevel!(dvs, ps, eqs, defs, network_sys::NeuronalNetworkSystem)
    
    return dvs, ps, eqs, defs, collect(neurons)
end

