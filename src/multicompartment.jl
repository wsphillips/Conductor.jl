
struct MultiCompartmentSystem <: AbstractCompartmentSystem
    iv::Num
    junctions::Set{<:AbstractJunction}
    compartments::Vector{CompartmentSystem}
    extensions::Set{ODESystem}
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
        return new(iv, compartments, extensions, eqs, systems, observed, defaults, name)
    end
end

const MultiCompartment = MultiCompartmentSystem

abstract type AbstractJunction end

struct Junction <: AbstractJunction
    parent::CompartmentSystem
    child::CompartmentSystem
    conductance::ConductanceSystem
    symmetric::Bool
end

Junction(x::Pair, cond; symmetric = true) = Junction(x.first, x.second, cond, symmetric)
get_conductance(x::Junction) = getfield(x, :conductance)
issymmetric(x::Junction) = getfield(x, :symmetric)

function MultiCompartment(junctions::Vector{<:AbstractJunction}; extensions = ODESystem[],
                          name = Base.gensym("MultiCompartment"))

    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    defaults = Dict()
    all_comp = Set{CompartmentSystem}()

    for jxn in junctions
        union!(all_comp, jxn.parent, jxn.child)
    end

    return MultiCompartment(t, junctions, all_comp, extensions, eqs, systems, observed,
                            defaults, name)
end

get_junctions(x::AbstractMultiCompartmentSystem) =  getfield(x, :junctions)
get_axial_conductance(x::AbstractCompartmentSystem) = getfield(x, :axial_conductance)
get_compartments(x::AbstractMultiCompartmentSystem) = getfield(x, :compartments)

MTK.get_states(x::AbstractMultiCompartmentSystem) = collect(build_toplevel(x)[1])
MTK.has_ps(x::AbstractMultiCompartmentSystem) = !isempty(build_toplevel(x)[2])
MTK.get_ps(x::AbstractMultiCompartmentSystem) = collect(build_toplevel(x)[2])

function MTK.get_eqs(x::AbstractMultiCompartmentSystem)
    empty!(getfield(x, :eqs))
    union!(getfield(x, :eqs), build_toplevel(x)[3])
end

MTK.get_defaults(x::AbstractMultiCompartmentSystem) = build_toplevel(x)[4]

function MTK.get_systems(x::AbstractMultiCompartmentSystem)
    empty!(getfield(x, :systems))
    union!(getfield(x, :systems), build_toplevel(x)[5], get_extensions(x))
end

function build_toplevel!(dvs, ps, eqs, defs, mcsys::MultiCompartmentSystem)

    junctions = get_junctions(mcsys)
    compartments = get_compartments(mcsys)
    jxn_types = Set{ConductanceSystem}() 

    forwards = Set{Equation}()

    for junction in junctions
        push!(jxn_types, get_conductance(junction))
    end
    
    # Reset subcompartment axial connections
    foreach(x -> empty!(get_axial_conductance(x)), compartments)
    
    for jxn in junctions
        axial = get_conductance(jxn) # maybe replicate? needs a toggle
        parent = jxn.parent
        child = jxn.child
        
        push!(get_axial_conductance(parent), axial)
        push!(forwards, child.Vₘ ~ getproperty(parent, nameof(axial)).Vₘ)

        if issymmetric(jxn)
            push!(get_axial_conductance(child), axial)
            push!(forwards, parent.Vₘ ~ getproperty(child, nameof(axial)).Vₘ)
        end           
    end

    union!(eqs, voltage_fwds)
    return dvs, ps, eqs, defs, collect(compartments)
end

#get_geometry(x::MultiCompartmentSystem)
#area(x::MultiCompartmentSystem)

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
    states, params, eqs, defs, compartments = build_toplevel(mcsys)

    all_systems = map(x -> convert(ODESystem, x), compartments)
    odesys = ODESystem(eqs, t, states, params; defaults = defs, name = nameof(mcsys))
    return compose(odesys, all_systems)
end
