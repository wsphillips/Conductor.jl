# Ion species
@enum IonSpecies::UInt128 begin
    NonIonic    = 1 << 0
    Sodium      = 1 << 1
    Potassium   = 1 << 2
    Chloride    = 1 << 3
    Calcium     = 1 << 4
end

const Ca = Calcium
const Na = Sodium
const K = Potassium
const Cl = Chloride
const Mixed = const Leak = NonIonic

const PERIODIC_SYMBOL = IdDict(Na => :Na, K  => :K, Cl => :Cl, Ca => :Ca, Leak => :l)

struct IonConcentration
    ion::IonSpecies
    loc::Location
end

const Concentration = IonConcentration

function IonConcentration(
    ion::IonSpecies,
    val = nothing;
    location::Location = Inside,
    dynamic::Bool = false,
    name::Symbol = PERIODIC_SYMBOL[ion]
)

    sym = Symbol(name,(loc == Inside ? "ᵢ" : "ₒ"))
    var = dynamic ? only(@variables $sym(t)) : only(@parameters $sym) 
    var = setmetadata(var,  IonConcentration, IonConcentration(ion, loc))
    if !isnothing(val)
        if val isa Molarity
            var = setmetadata(var, ConductorUnits, unit(val))
            raw_val = ustrip(Float64, val)
            var = setdefault(var, raw_val)
            return var
        else
            var = setdefault(var, val)
            return var
        end
    end
    return var
end

isconc(x) = hasmetadata(value(x), IonConcentration)
getconc(x) = isconc(x) ? getmetadata(value(x), IonConcentration) : nothing

struct IonCurrent
    ion::IonSpecies
    agg::Bool
end

function IonCurrent(
    ion::IonSpecies,
    val = nothing;
    aggregate::Bool = false,
    dynamic::Bool = false,
    name::Symbol = PERIODIC_SYMBOL[ion]
)
    sym = Symbol("I", name)
    var = dynamic ? only(@parameters $sym) : only(@variables $sym(t))
    setmetadata(var, IonCurrent, IonCurrent(ion, aggregate))
    if !isnothing(val)
        if val isa Current
            var = setmetadata(var, ConductorUnits, unit(val))
            raw_val = ustrip(Float64, val)
            var = setdefault(var, raw_val)
            return var
        else
            var = setdefault(var, val)
            return var
        end
    end
    return var
end

iscurrent(x) = hasmetadata(value(x), IonCurrent)
getcurrent(x) = iscurrent(x) ? getmetadata(value(x), IonCurrent) : nothing
getion(x::IonCurrent) = getfield(getcurrent(x), :ion)
isaggregate(x::IonCurrent) = getfield(getcurrent(x), :agg)

struct EquilibriumPotential
    ion::IonSpecies
end

const Equilibrium = EqulibirumPotential

function EquilibriumPotential(ion::IonSpecies, val; dynamic = false, name::Symbol = PERIODIC_SYMBOL[I])
    sym = Symbol("E", name)
    var = dynamic ? only(@variables $sym(t)) : only(@parameters $sym) 
    setmetadata(var, EquilibriumPotential, EquilibriumPotential(ion))
    if !isnothing(val)
        if val isa Voltage
            var = setmetadata(var, ConductorUnits, unit(val))
            raw_val = ustrip(Float64, val)
            var = setdefault(var, raw_val)
            return var
        else
            var = setdefault(var, val)
            return var
        end
    end
    return var
end

isreversal(x) = hasmetadata(value(x), EquilibriumPotential)
getreversal(x) = isreversal(x) ? getmetadata(value(x), EquilibriumPotential) : nothing
getion(x::EquilibriumPotential) = getfield(getreversal(x), :ion)

# Alternate constructor
function Equilibria(equil::Vector)
    out = Num[]
    for x in equil
        !(x.first <: IonSpecies) && throw("Equilibrium potential must be associated with an ion type.")
        if typeof(x.second) <: Tuple
            tup = x.second
            typeof(tup[2]) !== Symbol && throw("Second tuple argument for $(x.first) must be a symbol.")
            push!(out, Equilibrium(x.first, tup...))
        else
            push!(out, Equilibrium(x.first, x.second))
        end
    end
    return out
end
