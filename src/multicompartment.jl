
struct MultiCompartmentTopology
    g::SimpleDiGraph{Int}
    compartments::Vector{CompartmentSystem}
    conductances::Dict{Graphs.SimpleEdge{Int},ConductanceSystem}
end

function MultiCompartmentTopology(compartments::Vector{CompartmentSystem})
    g = SimpleDiGraph(length(compartments))
    conductances = Dict{Graphs.SimpleEdge{Int},ConductanceSystem}()
    return MultiCompartmentTopology(g, compartments, conductances)
end

nodes(topology::MultiCompartmentTopology) = getfield(topology, :compartments)
graph(topology::MultiCompartmentTopology) = getfield(topology, :g)

function find_compsys(compartment::AbstractCompartmentSystem, topology)
    return findfirst(isequal(compartment), vertices(topology))::Int 
end

function add_junction!(topology, trunk, branch, conductance::ConductanceSystem; symmetric = true) 
    src = find_compsys(trunk, topology)
    dst = find_compsys(branch, topology)
    add_edge!(graph(topology), src, dst)
    e = Graphs.SimpleEdge(src, dst)
    push!(topology.conductances, e => replicate(conductance))
    if symmetric
        add_edge!(topology.g, dst, src)
        e = Graphs.SimpleEdge(dst, src)
        push!(topology.conductances, e => replicate(conductance))
    end
    return nothing
end

function add_junction!(topology, x::Pair, conductance::ConductanceSystem; symmetric = true)
    add_junction!(topology, x.first, x.second, conductance, symmetric)
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
        foreach(x -> setparent!(x, mc), compartments)
        return mc
    end
end

const MultiCompartment = MultiCompartmentSystem
const NULL_AXIAL = Vector{Tuple{AbstractConductanceSystem,Num}}()

"""
$(TYPEDSIGNATURES)

Basic constructor for a `MultiCompartmentSystem`.
"""
function MultiCompartment(topology; extensions = ODESystem[],
                          name = Base.gensym("MultiCompartment"), defaults = Dict())
    
    compartments = topology.compartments
    
    # As a precaution, wipe any pre-existing axial currents
    for (i, comp) in enumerate(compartments)
        isempty(get_axial_conductance(comp)) && continue
        compartments[i] = CompartmentSystem(comp, axial_conductance = NULL_AXIAL)
    end

    observed = Equation[]
    eqs = Set{Equation}()
    
    for e in edges(topology.g)
        axial = topology.conductances[e]
        trunk = compartments[src(e)]
        branch = compartments[dst(e)]
        branchvm_alias = MembranePotential(nothing; name = Symbol(:V, nameof(branch))) 
        trunk = CompartmentSystem(trunk, axial_conductance = union(get_axial_conductance(trunk), [(axial, branchvm_alias)]))

        if hasproperty(getproperty(trunk, nameof(axial)), :Vₘ)
            push!(eqs, branch.Vₘ ~ getproperty(trunk, nameof(axial)).Vₘ)
        end
        
        push!(eqs, branch.Vₘ ~ getproperty(trunk, tosymbol(branchvm_alias, escape=false)))
        compartments[src(e)] = trunk
    end

    systems = union(extensions, compartments)

    return MultiCompartmentSystem(collect(eqs), t, [], [], observed, name, systems, defaults,
                                  topology, compartments, extensions) 
end

function MultiCompartmentSystem(x::MultiCompartmentSystem; topology = get_topology(x),
                                extensions = get_extensions(x), name = nameof(x),
                                defaults = get_defaults(x))
    MultiCompartmentSystem(topology;
                           extensions = extensions, name = name, defaults = defaults)
end

get_topology(x::MultiCompartmentSystem) = getfield(x, :topology)
get_compartments(x::MultiCompartmentTopology) = getfield(x, :compartments)
get_axial_conductance(x::CompartmentSystem) = getfield(x, :axial_conductance)
get_compartments(x::MultiCompartmentSystem) = getfield(x, :compartments)
get_compartments(x::CompartmentSystem) = [x]

function compartments(x::AbstractCompartmentSystem; namespace = true)
    compartment_systems = get_compartments(x) 
    if namespace && length(compartment_systems) > 1
        return [getproperty(x, nameof(n)) for n in compartment_systems]
    end
    return compartment_systems
end

hasparent(x::MultiCompartmentSystem) = false

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
    dvs = get_states(mcsys)
    ps  = get_ps(mcsys)
    eqs = get_eqs(mcsys)
    defs = get_defaults(mcsys)
    # why not get systems?
    systems = map(x -> convert(ODESystem, x), get_systems(mcsys))
    odesys = ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(mcsys), systems = systems)
    return odesys
end
