
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

function add_junction!(topology, trunk, branch, conductance::ConductanceSystem; symmetric = true) 
    src = findfirst(x -> isequal(x, trunk), topology.compartments)
    dst = findfirst(x -> isequal(x, branch), topology.compartments)
    (src === nothing || dst === nothing) && throw("junction compartments not found in topology")
    add_edge!(topology.g, src, dst)
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

get_axial_conductance(x::CompartmentSystem) = getfield(x, :axial_conductance)
get_compartments(x::MultiCompartmentSystem) = getfield(x, :compartments)
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
