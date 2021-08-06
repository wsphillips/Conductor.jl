using ExprTools, MacroTools, JuliaFormatter
using ExprTools: splitdef
using MacroTools: postwalk, prettify, prewalk, rmlines
using RuntimeGeneratedFunctions: @RuntimeGeneratedFunction

isivsuffix(x,sys) = endswith(string(x),"($(independent_variable(sys)))")

function chopiv(x,sys)
    suflen = length(string(independent_variable(sys))) + 2
    return Symbol(chop(string(x), tail = suflen))
end

function rgf_lambda_expr(x)
    args = first(ExprTools.parameters(typeof(x)))
    return Expr(:macrocall, Symbol("@RuntimeGeneratedFunction"), Expr(:function, Expr(:tuple, args...), x.body))
end

function clean_expr!(ex_org, sys; rgf=false)
    ex = prewalk(rmlines, ex_org)
    def = splitdef(ex)
    for (old, new) in zip(def[:args][1:end-1], [:du, :u, :p])
        ex = postwalk(x -> x isa Symbol && x == old ? new : x, ex)
    end
    if rgf
        ex = postwalk(x -> x isa RuntimeGeneratedFunction ? rgf_lambda_expr(x) : x, ex)
        ex = prewalk(rmlines, ex)
    end
    ex = postwalk(x -> (x isa Symbol && isivsuffix(x,sys)) ? chopiv(x,sys) : x, ex)
    return ex
end


