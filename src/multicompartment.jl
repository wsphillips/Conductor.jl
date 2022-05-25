
struct MultiCompartmentSystem <: AbstractCompartmentSystem
    iv::Num
    junctions::Set{Junction}
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

struct ScaledJunction <: AbstractJunction
    compartments::Vector{AbstractCompartmentSystem}
    scale_factor::Num
end

function MultiCompartment(junctions::Vector{<:AbstractJunction}; extensions = ODESystem[],
                          name = Base.gensym("MultiCompartment"))

    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    defaults = Dict()
    all_comp = Set{CompartmentSystem}()

    for jxn in junctions
        union!(all_comp, jxn.compartments)
    end

    return MultiCompartment(t, junctions, all_comp, extensions, eqs, systems, observed,
                            defaults, name)
end

function build_toplevel!(mcsys::MultiCompartmentSystem)
end

get_geometry(x::MultiCompartmentSystem)
area(x::MultiCompartmentSystem)

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
end
