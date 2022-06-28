abstract type AbstractJunction end

"""
$(TYPEDEF)

A connection (edge) between two morphologically contiguous compartments.

For example, a `Junction` between a somatic `CompartmentSystem` and a dendritic `CompartmentSystem`.
"""
struct Junction <: AbstractJunction
    trunk::CompartmentSystem
    branch::CompartmentSystem
    conductance::ConductanceSystem
    symmetric::Bool
end

"""
    $(TYPEDSIGNATURES)

Basic constructor for a `Junction`. 

"""
function Junction(x::Pair, conductance::ConductanceSystem; symmetric = true)
    return Junction(x.first, x.second, conductance, symmetric)
end

get_conductance(x::Junction) = getfield(x, :conductance)
issymmetric(x::Junction) = getfield(x, :symmetric)

"""
$(TYPEDEF)

A neuron with 2+ morphologically connected compartments.

$(TYPEDFIELDS)
"""
struct MultiCompartmentSystem <: AbstractCompartmentSystem
    "Independent variabe. Defaults to time, ``t``."
    iv::Num
    "`Junction` (edges) between subcompartments of the neuron."
    junctions::Vector{AbstractJunction}
    "Individual subcompartments of the neuron."
    compartments::Vector{CompartmentSystem}
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    eqs::Vector{Equation}
    systems::Vector{AbstractTimeDependentSystem}
    observed::Vector{Equation}
    defaults::Dict
    name::Symbol
    function MultiCompartmentSystem(iv, junctions, compartments, extensions, eqs, systems,
                                    observed, defaults, name; checks = false) 
        if checks
            #placeholder
        end
        mc = new(iv, junctions, compartments, extensions, eqs, systems, observed, defaults,
                 name)
        foreach(x -> setparent!(x, mc), compartments)
        return mc
    end
end

const MultiCompartment = MultiCompartmentSystem

"""
$(TYPEDSIGNATURES)

Basic constructor for a `MultiCompartmentSystem`.
"""
function MultiCompartment(junctions::Vector{<:AbstractJunction}; extensions = ODESystem[],
                          name = Base.gensym("MultiCompartment"))

    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    defaults = Dict()
    all_comp = Set{CompartmentSystem}()

    for jxn in junctions
        push!(all_comp, jxn.branch)
        push!(all_comp, jxn.trunk)
    end
    
    return MultiCompartment(t, junctions, collect(all_comp), extensions, eqs, systems,
                            observed, defaults, name)
end

get_junctions(x::MultiCompartmentSystem) =  getfield(x, :junctions)
get_axial_conductance(x::AbstractCompartmentSystem) = getfield(x, :axial_conductance)
get_compartments(x::MultiCompartmentSystem) = getfield(x, :compartments)
hasparent(x::MultiCompartmentSystem) = false

function get_systems(x::MultiCompartmentSystem; rebuild = false)
    empty!(getfield(x, :systems))
    union!(getfield(x, :systems), getfield(x, :compartments), getfield(x, :extensions))
    return getfield(x, :systems)
end

function build_toplevel!(dvs, ps, eqs, defs, mcsys::MultiCompartmentSystem)

    junctions = get_junctions(mcsys)
    compartments = get_compartments(mcsys)
    forwards = Set{Equation}()

    # Reset subcompartment axial connections
    foreach(x -> empty!(get_axial_conductance(x)), compartments)
    
    for jxn in junctions
        axial = get_conductance(jxn)
        trunk = jxn.trunk
        branch = jxn.branch
        branchvm = MembranePotential(; name = Symbol(:V, nameof(branch))) 
        push!(get_axial_conductance(trunk), (axial, branchvm))
        #TODO: resolve arbitrary states generically, not just hardcoded Vₘ
        if hasproperty(getproperty(trunk, nameof(axial)), :Vₘ)
            push!(forwards, branch.Vₘ ~ getproperty(trunk, nameof(axial)).Vₘ)
        end
        push!(forwards, branch.Vₘ ~ getproperty(trunk, tosymbol(branchvm, escape=false)))
        if issymmetric(jxn)
            trunkvm = MembranePotential(; name = Symbol(:V, nameof(trunk)))
            push!(get_axial_conductance(branch), (axial, trunkvm))
            #TODO: resolve arbitrary states generically, not just hardcoded Vₘ
            if hasproperty(getproperty(branch, nameof(axial)), :Vₘ)
                push!(forwards, trunk.Vₘ ~ getproperty(branch, nameof(axial)).Vₘ)
            end
            push!(forwards, trunk.Vₘ ~ getproperty(branch, tosymbol(trunkvm, escape=false)))
        end
    end

    union!(eqs, forwards)
    return dvs, ps, eqs, defs, collect(compartments)
end

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
    states, params, eqs, defs, compartments = build_toplevel(mcsys)
    all_systems = map(x -> convert(ODESystem, x), compartments)
    odesys = ODESystem(eqs, t, states, params; defaults = defs, name = nameof(mcsys), systems = all_systems)
    return odesys
end
