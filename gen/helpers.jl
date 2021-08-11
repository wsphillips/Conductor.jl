using ExprTools, MacroTools, JuliaFormatter
using ExprTools: splitdef
using MacroTools: postwalk, prettify, prewalk, rmlines
using RuntimeGeneratedFunctions

isivsuffix(x,sys) = endswith(string(x),"($(independent_variable(sys)))")

function chopiv(x,sys)
    suflen = length(string(independent_variable(sys))) + 2
    return Symbol(chop(string(x), tail = suflen))
end

function rgf_lambda_expr(x)
    argnames,cache_tag,context_tag,id = ExprTools.parameters(typeof(x))
    return :( @RuntimeGeneratedFunction($(Expr(:quote,:(($(argnames...),)->$(x.body)) ))))
end

function clean_expr(ex_org, sys; rgf=false, pretty=false, name=nothing)
    def = splitdef(ex_org)
    body = prewalk(rmlines, def[:body])
    for (old, new) in zip(def[:args][1:3], [:du, :u, :p])
        body = postwalk(x -> x isa Symbol && x == old ? new : x, body)
    end
    if rgf
        body = postwalk(x -> x isa RuntimeGeneratedFunction ? rgf_lambda_expr(x) : x, body)
        body = prewalk(rmlines, body)
    end
    body = postwalk(x -> (x isa Symbol && isivsuffix(x,sys)) ? chopiv(x,sys) : x, body)
    
    def[:args][1:3] = [:du, :u, :p]
    def[:body] = body
    if !isnothing(name)
        def[:name] = name
    end
    return pretty ? prettify(ExprTools.combinedef(def)) : ExprTools.combinedef(def)
end

display_expr(ex) = println(format_text(string(ex)))

function write_expr(filename, ex)
    open(filename, "w") do io
        write(io, format_text(string(ex)))
    end
end
