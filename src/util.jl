
using IfElse
import Symbolics: unwrap, symtype, getindex_posthook

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

const ExprValues = Union{Expr,Symbol,Number}  # For use in macros
isfunction(ex::ExprValues) = try return eval(ex) isa Function catch; return false end

function extract_symbols(ex::ExprValues, out::Vector{Symbol}=[])
    if ~isfunction(ex) && isa(ex, Symbol)
        union!(out, [ex])
    end
    return ex
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

