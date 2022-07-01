
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

namegen(name) = Symbol(filter(x -> x !== '#', String(Base.gensym(name))))

function replicate(x::Union{AbstractCompartmentSystem,AbstractConductanceSystem})
    rootname = ModelingToolkit.getname(x)
    new = deepcopy(x)
    return ModelingToolkit.rename(new, namegen(rootname))
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

function NetworkTopology(neurons::Vector{<:AbstractCompartmentSystem},
                         synaptic_models::Vector{<:AbstractConductanceSystem})
    n = length(neurons)
    multigraph = Dict([x => sparse(Int64[], Int64[], Num[], n, n) for x in synaptic_models])
    return NetworkTopology(multigraph, neurons)
end

function NetworkTopology(g::SimpleDiGraph, neurons::Vector{<:AbstractCompartmentSystem},
        synaptic_model::AbstractConductanceSystem, default_weight = get_gbar(synaptic_model))
    compartments = CompartmentSystem[]

    for neuron in neurons
          
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

  #  for synaptic_class in keys(multigraph)
  #      
  #      # fetch the adjacency matrix
  #      g = multigraph[synaptic_class]
  #      rows = rowvals(g)
  #      vals = nonzeros(g)
  #      for i in axes(g, 2)
  #          post_compartment =             
  #          pre_compartments = 

    preneurons  = Set()
    postneurons = Set()
    neurons = Set()
    synaptic_channel_classes = Set()

    # Bin the fields
    for synapse in synapses
        push!(preneurons, presynaptic(synapse))
        push!(postneurons, postsynaptic(synapse))
        push!(synaptic_channel_classes, class(synapse))
    end
    
    all_compartments = union(preneurons, postneurons)
    
    neurons = deepcopy.(neurons)
    postneurons = deepcopy.(postneurons)
    # For each post synaptic neuron
    for pn in postneurons
        incoming_edges = filter(x -> isequal(postsynaptic(x), pn), synapses)
        for cl in synaptic_channel_classes
            # filter to incoming synapses of the same type per postsynaptic compartment
            inc_edges_of_cl = filter(x -> isequal(class(x), cl), incoming_edges)
            cl_copy = deepcopy(cl) 
            for syn in inc_edges_of_cl
                if !isaggregate(cl_copy)
                    cl_copy = replicate(class(syn))
                end
                subscribe!(cl_copy, presynaptic(syn))
                push!(get_synaptic_reversals(pn), reversal(syn))
                if !isaggregate(cl_copy)
                    push!(get_synapses(pn), cl_copy)
                    post_v = getproperty(pn, nameof(cl_copy)).Vₓ
                    pre_v = presynaptic(syn).Vₘ
                    push!(voltage_fwds, pre_v ~ post_v)
                    push!(defs, post_v => pre_v)
                end
            end
            if !isempty(subscriptions(cl_copy)) && isaggregate(cl_copy)
                push!(get_synapses(pn), cl_copy)

                if length(inc_edges_of_cl) > 1
                    vxs = filter(x -> isvoltage(x) && isextrinsic(x), states(getproperty(pn, nameof(cl_copy))))
                    for (vx, pre) in zip(vxs, presynaptic.(inc_edges_of_cl))
                        post_v = renamespace(getproperty(pn, nameof(cl_copy)), vx)
                        pre_v = pre.Vₘ
                        push!(voltage_fwds, pre_v ~ post_v)
                        push!(defs, post_v => pre_v)
                    end
                else
                    post_v = getproperty(pn, nameof(cl_copy)).Vₓ
                    pre_v = only(presynaptic.(inc_edges_of_cl)).Vₘ
                    push!(voltage_fwds, pre_v ~ post_v)
                    push!(defs, post_v => pre_v)
                end
            end
        end
    end

    union!(eqs, voltage_fwds)



    return NeuronalNetworkSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                                   topology, reversal_map; checks = false)
end

get_extensions(x::AbstractNeuronalNetworkSystem) = getfield(x, :extensions)
get_synapses(x::AbstractNeuronalNetworkSystem) = getfield(x, :synapses)

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

