
using IfElse

function heaviside(x)
    IfElse.ifelse(x > zero(x), one(x), zero(x))
end
