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

struct MultiCompartmentSystem <: AbstractCompartmentSystem
    iv::Num
    junctions::Vector{AbstractJunction}
    compartments::Vector{CompartmentSystem}
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

function MultiCompartment(junctions::Vector{<:AbstractJunction}; extensions = ODESystem[],
                          name = Base.gensym("MultiCompartment"))

    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    defaults = Dict()
    all_comp = Set{CompartmentSystem}()

    for jxn in junctions
        push!(all_comp, jxn.child)
        push!(all_comp, jxn.parent)
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
    jxn_types = Set{ConductanceSystem}() 

    forwards = Set{Equation}()

    for junction in junctions
        push!(jxn_types, get_conductance(junction))
    end
    
    # Reset subcompartment axial connections
    foreach(x -> empty!(get_axial_conductance(x)), compartments)
    
    for jxn in junctions
        axial = get_conductance(jxn)
        parent = jxn.parent
        child = jxn.child
        childvm = MembranePotential(; name = Symbol(:V, nameof(child))) 
        push!(get_axial_conductance(parent), (axial, childvm))
        #TODO: resolve arbitrary states generically
        if hasproperty(getproperty(parent, nameof(axial)), :Vₘ)
            push!(forwards, child.Vₘ ~ getproperty(parent, nameof(axial)).Vₘ)
        end
        push!(forwards, child.Vₘ ~ getproperty(parent, tosymbol(childvm, escape=false)))
        if issymmetric(jxn)
            parentvm = MembranePotential(; name = Symbol(:V, nameof(parent)))
            push!(get_axial_conductance(child), (axial, parentvm))
            #TODO: resolve arbitrary states generically
            if hasproperty(getproperty(child, nameof(axial)), :Vₘ)
                push!(forwards, parent.Vₘ ~ getproperty(child, nameof(axial)).Vₘ)
            end
            push!(forwards, parent.Vₘ ~ getproperty(child, tosymbol(parentvm, escape=false)))
        end
    end

    union!(eqs, forwards)
    return dvs, ps, eqs, defs, collect(compartments)
end

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
    states, params, eqs, defs, compartments = build_toplevel(mcsys)
    all_systems = map(x -> convert(ODESystem, x), compartments)
    odesys = ODESystem(eqs, t, states, params; defaults = defs, name = nameof(mcsys))
    return compose(odesys, all_systems)
end
