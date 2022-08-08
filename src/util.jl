
using IfElse
import Symbolics: unwrap, symtype, getindex_posthook

namegen(name) = Symbol(filter(x -> x !== '#', String(Base.gensym(name))))

function replicate(x::Union{AbstractCompartmentSystem,AbstractConductanceSystem})
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

