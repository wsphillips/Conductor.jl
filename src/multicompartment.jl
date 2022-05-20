
struct MultiCompartmentSystem <: AbstractCompartmentSystem
    iv::Num
    junctions::Set{Junction}
    compartments::Set{CompartmentSystem}
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

struct Junction
    compartments::Vector{AbstractCompartmentSystem}
end

function MultiCompartment(jxns::Vector{Junction}; extensions = ODESystem[],
                          name = Base.gensym("MultiCompartment"))

    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    defaults = Dict()



    return MultiCompartment(t, jxns, compartments, extensions, eqs, systems, observed,
                            defaults, name)
end

function build_toplevel!(mcsys::MultiCompartmentSystem)
end

get_geometry(x::MultiCompartmentSystem)
area(x::MultiCompartmentSystem)

function Base.convert(::Type{ODESystem}, mcsys::MultiCompartmentSystem)
end
