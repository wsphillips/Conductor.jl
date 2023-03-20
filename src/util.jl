# basic fields needed for AbstractSystem interface in MTK
struct GeneratedCollections
    eqs::Vector{Equation}
    dvs::Set{Num}
    ps::Set{Num}
    systems::Vector{AbstractTimeDependentSystem}
    observed::Vector{Equation}
    defs::Dict
end

function GeneratedCollections(; eqs = Equation[], systems = AbstractTimeDependentSystem[],
                              observed = Equation[], dvs = Set{Num}(), ps = Set{Num}(),
                              defs = Dict())
    return GeneratedCollections(eqs, dvs, ps, systems, observed, defs)
end

function copy_collections(sys::AbstractSystem) 
    GeneratedCollections(eqs = deepcopy(get_eqs(sys)),
                         systems = deepcopy(get_systems(sys)),
                         observed = deepcopy(get_observed(sys)),
                         dvs = Set(deepcopy(get_states(sys))),
                         ps = Set(deepcopy(get_ps(sys))),
                         defs = deepcopy(get_defaults(sys)))
end

import Symbolics: unwrap, symtype, getindex_posthook

indexof(sym, syms) = findfirst(isequal(sym), syms)

function indexmap(syms, ref)
    idxs = Int64[]
    for sym in syms
        idx = indexof(sym, ref)
        isnothing(idx) && continue # drop symbols that dont exist
        push!(idxs, idx)
    end
    return idxs
end

namegen(name) = Symbol(filter(x -> x !== '#', String(Base.gensym(name))))

function replicate(x::Union{AbstractCompartmentSystem, AbstractConductanceSystem})
    rootname = ModelingToolkit.getname(x)
    new = deepcopy(x)
    return ModelingToolkit.rename(new, namegen(rootname))
end

function build_toplevel(system)
    dvs = Set{Num}()
    ps = Set{Num}()
    eqs = Equation[]
    defs = Dict()
    build_toplevel!(dvs, ps, eqs, defs, system)
end

function time_span(tspan::Tuple{Time,Time})
    ustrip(Float64, ms, tspan[1]), ustrip(Float64, ms, tspan[2])
end

time_span(tspan::Tuple{Float64,Float64}) = tspan
time_span(tspan::Time) = zero(Float64), ustrip(Float64, ms, tspan)
time_span(tspan::Real) = zero(Float64), Float64(tspan)

heaviside(x) = ifelse(x > zero(x), one(x), zero(x))
@register_symbolic heaviside(x)
ModelingToolkit.get_unit(op::typeof(heaviside), args) = ms^-1

# Hijacked and modified from Symbolics.jl
function set_symarray_metadata(x, ctx, val)
    if symtype(x) <: AbstractArray
        if val isa AbstractArray
            getindex_posthook(x) do r, x, i...
                set_symarray_metadata(r, ctx, val[i...])
            end
        else
            getindex_posthook(x) do r, x, i...
                set_symarray_metadata(r, ctx, val)
            end
        end
    else
        setmetadata(x, ctx, val)
    end
end

function settunable(sys, ps)
    new_ps = setmetadata.(ps, ModelingToolkit.VariableTunable, true) 
    @set! sys.ps = union(new_ps, parameters(sys))
    return sys
end

function setbounds(sys, p_dists)
    new_ps = parameters(sys)
    for pair in p_dists
        i = findfirst(isequal(first(pair)), new_ps)
        i == nothing && throw("Parameter $(first(pair)) not found.")
        new_ps[i] = setmetadata.(new_ps, ModelingToolkit.VariableBounds, second(pair))
    end
    @set! sys.ps = union(new_ps, parameters(sys))
    return sys
end

function setdistributions(sys, ps) end

