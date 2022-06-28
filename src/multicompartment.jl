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
    "`Junction` (edges) between subcompartments of the neuron."
    junctions::Vector{AbstractJunction}
    "Individual subcompartments of the neuron."
    compartments::Vector{CompartmentSystem}
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    function MultiCompartmentSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                                    junctions, compartments, extensions; checks = false) 
        if checks
            #placeholder
        end
        mc = new(eqs, iv, states, ps, observed, name, systems, defaults,
                 junctions, compartments, extensions)
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
        name = Base.gensym("MultiCompartment"), defaults = Dict())

    compartments = Set{CompartmentSystem}()

    for jxn in junctions
        push!(compartments, jxn.branch)
        push!(compartments, jxn.trunk)
    end
    
    return MultiCompartment(junctions, collect(compartments), extensions, eqs, systems,
                            observed, defaults, name)
end

get_junctions(x::MultiCompartmentSystem) =  getfield(x, :junctions)
get_axial_conductance(x::AbstractCompartmentSystem) = getfield(x, :axial_conductance)
get_compartments(x::MultiCompartmentSystem) = getfield(x, :compartments)
hasparent(x::MultiCompartmentSystem) = false

function MultiCompartment(junctions, compartments, extensions, defaults, name)

    systems = AbstractTimeDependentSystem[]
    observed = Equation[]

    eqs = Set{Equation}()

    # MUTATION Reset subcompartment axial connections
    foreach(x -> empty!(get_axial_conductance(x)), compartments)
    
    for jxn in junctions
        axial = get_conductance(jxn)
        trunk = jxn.trunk
        branch = jxn.branch
        branchvm = MembranePotential(; name = Symbol(:V, nameof(branch))) 
        # When this gets updated we should rebuild the compartment system
        # perhaps overload `setproperties` for SetField.jl
        push!(get_axial_conductance(trunk), (axial, branchvm))

        #TODO: resolve arbitrary states generically, not just hardcoded Vₘ
        # connect (branch compartment -> vm) to (trunk comp -> axial conductance -> vm)
        if hasproperty(getproperty(trunk, nameof(axial)), :Vₘ)
            push!(eqs, branch.Vₘ ~ getproperty(trunk, nameof(axial)).Vₘ)
        end

        push!(eqs, branch.Vₘ ~ getproperty(trunk, tosymbol(branchvm, escape=false)))

        if issymmetric(jxn)
            trunkvm = MembranePotential(; name = Symbol(:V, nameof(trunk)))

            ### mutation
            push!(get_axial_conductance(branch), (axial, trunkvm))
            ### mutation
            #
            #TODO: resolve arbitrary states generically, not just hardcoded Vₘ
            if hasproperty(getproperty(branch, nameof(axial)), :Vₘ)
                push!(eqs, trunk.Vₘ ~ getproperty(branch, nameof(axial)).Vₘ)
            end
            push!(eqs, trunk.Vₘ ~ getproperty(branch, tosymbol(trunkvm, escape=false)))
        end
    end

    return dvs, ps, eqs, defs, collect(compartments)
end

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
    dvs = get_states(mcsys)
    ps  = get_ps(mcsys)
    eqs = get_eqs(mcsys)
    defs = get_defaults(mcsys)
    compartments = get_compartments(mcsys)
    all_systems = map(x -> convert(ODESystem, x), compartments)
    odesys = ODESystem(eqs, t, states, params; defaults = defs, name = nameof(mcsys), systems = all_systems)
    return odesys
end
