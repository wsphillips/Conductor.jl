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
                              observed = Equation[], dvs = Set{Num}(), ps  = Set{Num}(),
                              defs = Dict())
    return GeneratedCollections(eqs, dvs, ps, systems, observed, defs)
end

using IfElse
import Symbolics: unwrap, symtype, getindex_posthook

namegen(name) = Symbol(filter(x -> x !== '#', String(Base.gensym(name))))

function replicate(x::Union{AbstractCompartmentSystem,AbstractConductanceSystem})
    rootname = ModelingToolkit.getname(x)
    new = deepcopy(x)
    return ModelingToolkit.rename(new, namegen(rootname))
end

function genvar(name; iv = nothing, default = nothing)
    x = name
    if isnothing(default)
        return only(isnothing(iv) ? @parameters($x) : @variables($x(iv)))
    else
        return only(isnothing(iv) ? @parameters($x=default) :  @variables($x(iv)=default))
    end
end

function build_toplevel(system)
    dvs = Set{Num}()
    ps = Set{Num}()
    eqs = Equation[]
    defs = Dict()
    build_toplevel!(dvs, ps, eqs, defs, system)
end

function heaviside(x)
    IfElse.ifelse(x > zero(x), one(x), zero(x))
end

# Hijacked and modified from Symbolics.jl
function set_symarray_metadata(x, ctx, val)
    if symtype(x) <: AbstractArray
        if val isa AbstractArray
            getindex_posthook(x) do r,x,i...
                set_symarray_metadata(r, ctx, val[i...])
            end
        else
            getindex_posthook(x) do r,x,i...
                set_symarray_metadata(r, ctx, val)
            end
        end
    else
        setmetadata(x, ctx, val)
    end
end

#=
const ExprValues = Union{Expr,Symbol,Number}  # For use in macros
isfunction(ex::ExprValues) = try return eval(ex) isa Function catch; return false end

function extract_symbols(ex::ExprValues, out::Vector{Symbol}=[])
    if ~isfunction(ex) && isa(ex, Symbol)
        union!(out, [ex])
    end
    return ex
end
=#

