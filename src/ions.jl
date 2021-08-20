# Metadata IDs
abstract type ConductorCurrentCtx end
abstract type ConductorEquilibriumCtx end
abstract type ConductorConcentrationCtx end
abstract type ConductorAggregatorCtx end

# Ion species
abstract type Ion end
abstract type Cation <: Ion end
abstract type Anion <: Ion end
struct Calcium  <: Ion end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Chloride <: Ion end

const Ca = Calcium
const Na = Sodium
const K = Potassium
const Cl = Chloride
const Mixed = Ion # non-specific ion
const Leak = Mixed

const PERIODIC_SYMBOL = IdDict(Na => :Na, K  => :K, Cl => :Cl, Ca => :Ca, Leak => :l)

# Concentrations of ions
abstract type AbstractConcentration end

struct IonConcentration{I<:Ion, L<:Location, V<:Union{Nothing, Num,Symbolic,Molarity}} <:AbstractConcentration
    ion::Type{I}
    val::V
    loc::L
end

# FIXME: handle default values better
function Concentration(::Type{I}, val = nothing, loc::Location = Inside, name::Symbol = PERIODIC_SYMBOL[I]) where {I <: Ion}
    sym = Symbol(name,(loc == Inside ? "ᵢ" : "ₒ"))
    # FIXME: Not necessarily a parameter when used as a primitive...but we should support
    # this. Use a flag?
    var = #=val isa Molarity ? only(@parameters $sym) :=# only(@variables $sym(t))
    return setmetadata(var,  ConductorConcentrationCtx, IonConcentration(I, val, loc))
end

isconcentration(x::Symbolic) = hasmetadata(x, ConductorConcentrationCtx)
isconcentration(x::Num) = isconcentration(value(x))
getconcentration(x::Symbolic) = isconcentration(x) ? getmetadata(x, ConductorConcentrationCtx) : nothing
getconcentration(x::Num) = getconcentration(value(x))

# Currents
struct MembraneCurrent{I<:Ion,V<:Union{Nothing,Num,Symbolic,Current}}
    ion::Type{I}
    val::V
end

# TODO: add aggregator as field of Membrane current struct
function MembraneCurrent{I}(val = nothing; name::Symbol = PERIODIC_SYMBOL[I], aggregate::Bool = false) where {I <: Ion}
    sym = Symbol("I", name)
    var = val isa Current ? only(@parameters $sym) : only(@variables $sym(t))
    var = setmetadata(var, ConductorCurrentCtx, MembraneCurrent(I, val))
    return setmetadata(var, ConductorAggregatorCtx, aggregate)
end

ismembranecurrent(x::Symbolic) = hasmetadata(x, ConductorCurrentCtx)
ismembranecurrent(x::Num) = ismembranecurrent(ModelingToolkit.value(x))
getmembranecurrent(x::Union{Num, Symbolic}) = ismembranecurrent(x) ? getmetadata(x, ConductorCurrentCtx) : nothing
iontype(x::Union{Num, Symbolic}) = getmembranecurrent(x).ion
isaggregator(x::Union{Num, Symbolic})  = getmetadata(x, ConductorAggregatorCtx)

# Equilibrium potential implicitly defines an ionic gradient
abstract type AbstractIonGradient end

struct EquilibriumPotential{I<:Ion,V<:Union{Num,Symbolic,Voltage}} <: AbstractIonGradient
    ion::Type{I}
    val::V
end

const Equilibrium{I} = EquilibriumPotential{I}

function EquilibriumPotential{I}(val, name::Symbol = PERIODIC_SYMBOL[I]) where {I <: Ion}
    sym = Symbol("E", name)
    var = val isa Voltage ? only(@parameters $sym) : only(@variables $sym(t))
    return setmetadata(var, ConductorEquilibriumCtx, EquilibriumPotential(I, val))
end

# Alternate constructor
function Equilibria(equil::Vector)
    out = Num[]
    for x in equil
        !(x.first <: Ion) && throw("Equilibrium potential must be associated with an ion type.")
        if typeof(x.second) <: Tuple
            tup = x.second
            typeof(tup[2]) !== Symbol && throw("Second tuple argument for $(x.first) must be a symbol.")
            push!(out, Equilibrium{x.first}(tup...))
        else
            push!(out, Equilibrium{x.first}(x.second))
        end
    end
    return out
end
