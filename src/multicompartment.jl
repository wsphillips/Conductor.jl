
struct MultiCompartmentTopology
    g::SimpleDiGraph{Int}
    compartments::Vector{CompartmentSystem}
    conductances::Dict{Graphs.SimpleEdge{Int}, ConductanceSystem}
end

function MultiCompartmentTopology(compartments::Vector{CompartmentSystem{T}}) where {T <: AbstractDynamics}
    g = SimpleDiGraph(length(compartments))
    conductances = Dict{Graphs.SimpleEdge{Int}, ConductanceSystem}()
    return MultiCompartmentTopology(g, compartments, conductances)
end

vertices(topology::MultiCompartmentTopology) = getfield(topology, :compartments)
graph(topology::MultiCompartmentTopology) = getfield(topology, :g)

function find_compsys(compartment::AbstractCompartmentSystem, topology)
    return findfirst(isequal(nameof(compartment)), nameof.(vertices(topology)))::Int
end

function add_junction!(topology, trunk, branch, conductance::ConductanceSystem)
    src = find_compsys(trunk, topology)
    dst = find_compsys(branch, topology)

    add_edge!(graph(topology), src, dst)
    e = Graphs.SimpleEdge(src, dst)
    push!(topology.conductances, e => replicate(conductance))

    add_edge!(topology.g, dst, src)
    e = Graphs.SimpleEdge(dst, src)
    push!(topology.conductances, e => replicate(conductance))
    return nothing
end

function add_junction!(topology, trunk, branch, conductances::NTuple{2, ConductanceSystem})
    src = find_compsys(trunk, topology)
    dst = find_compsys(branch, topology)

    add_edge!(graph(topology), src, dst)
    e = Graphs.SimpleEdge(src, dst)
    push!(topology.conductances, e => replicate(conductances[1]))

    add_edge!(topology.g, dst, src)
    e = Graphs.SimpleEdge(dst, src)
    push!(topology.conductances, e => replicate(conductances[2]))
    return nothing
end

function add_junction!(topology, x::Pair, conductances::NTuple{2, ConductanceSystem})
    add_junction!(topology, x.first, x.second, conductances)
end

"""
$(TYPEDEF)

A neuron with 2+ morphologically connected compartments.

$(TYPEDFIELDS)
"""
struct MultiCompartmentSystem <: AbstractCompartmentSystem
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
    topology::MultiCompartmentTopology
    "Individual subcompartments of the neuron."
    compartments::Vector{CompartmentSystem}
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    function MultiCompartmentSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                                    topology, compartments, extensions; checks = false)
        if checks
            #placeholder
        end
        mc = new(eqs, iv, states, ps, observed, name, systems, defaults,
                 topology, compartments, extensions)
        return mc
    end
end

const MultiCompartment = MultiCompartmentSystem
const NULL_AXIAL = Vector{Tuple{AbstractConductanceSystem, Num}}()

"""
$(TYPEDSIGNATURES)

Basic constructor for a `MultiCompartmentSystem`.
"""
function MultiCompartment(topology::MultiCompartmentTopology; extensions = ODESystem[],
                          name = Base.gensym("MultiCompartment"), defaults = Dict())
    compartments = topology.compartments

    # As a precaution, wipe any pre-existing axial currents from compartments
    for (i, comp) in enumerate(compartments)
        isempty(conductances(get_arbor(comp))) && continue
        compartments[i] = SciMLBase.remake(comp, arbor = Arborization())
    end

    observed = Equation[]
    eqs = Set{Equation}()

    for e in edges(topology.g)
        axial = topology.conductances[e]
        trunk = compartments[src(e)]
        branch = compartments[dst(e)]

        trunk_arbor = get_arbor(trunk)
        push!(trunk_arbor.children,
              Junction(axial, Num(renamespace(branch, get_voltage(branch)))))

        trunk = SciMLBase.remake(trunk, arbor = trunk_arbor)
        compartments[src(e)] = trunk
    end

    systems = union(extensions, compartments)

    return MultiCompartmentSystem(collect(eqs), t, [], [], observed, name, systems,
                                  defaults,
                                  topology, compartments, extensions)
end

function SciMLBase.remake(x::MultiCompartmentSystem; topology = get_topology(x),
                          extensions = get_extensions(x), name = nameof(x),
                          defaults = get_defaults(x))
    MultiCompartmentSystem(topology;
                           extensions = extensions, name = name, defaults = defaults)
end

get_topology(x::MultiCompartmentSystem) = getfield(x, :topology)
get_compartments(x::MultiCompartmentTopology) = getfield(x, :compartments)
get_compartments(x::MultiCompartmentSystem) = getfield(x, :compartments)
get_compartments(x::CompartmentSystem) = [x]

function compartments(x::AbstractCompartmentSystem; namespace = true)
    compartment_systems = get_compartments(x)
    if namespace && length(compartment_systems) > 1
        return [getproperty(x, nameof(n)) for n in compartment_systems]
    end
    return compartment_systems
end

Base.eltype(::MultiCompartmentSystem) = CompartmentSystem
Base.length(M::MultiCompartmentSystem) = length(get_compartments(M))

function Base.iterate(M::MultiCompartmentSystem, state = 1)
    state > length(M) && return nothing
    return (get_compartments(M)[state], state + 1)
end

function Base.iterate(rM::Iterators.Reverse{MultiCompartmentSystem}, state = length(rM.itr))
    state < 1 && return nothing
    mc = rM.itr
    return (get_compartments(mc)[state], state - 1)
end

function Base.getindex(M::MultiCompartmentSystem, i; namespace = true)
    if namespace
        return getproperty(M, nameof(get_compartments(M)[i]))
    else
        return get_compartments(M)[i]
    end
end

Base.firstindex(M::MultiCompartmentSystem) = 1
Base.lastindex(M::MultiCompartmentSystem) = length(M)

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
    dvs = get_states(mcsys)
    ps = get_ps(mcsys)
    eqs = get_eqs(mcsys)
    defs = get_defaults(mcsys)
    systems = map(x -> convert(ODESystem, x), get_systems(mcsys))
    odesys = ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(mcsys),
                       systems = systems, checks = CheckComponents)
    return odesys
end
