const ℱ = Unitful.q*Unitful.Na # Faraday's constant
const t = let name = :t; only(@variables $name) end
const D = Differential(t)

# Properties
@enum PrimitiveSource Intrinsic Extrinsic
@enum PrimitiveLocation Outside Inside

struct MembranePotential
    function MembranePotential(V0 = -60mV; dynamic = true, source::PrimitiveSource = Intrinsic,
                               n::Integer = 1, name::Symbol = :Vₘ)
        if isnothing(V0)
            if n == one(n) 
                ret = dynamic ? only(@variables $name(t)) :
                                only(@parameters $name)
            elseif n > one(n)
                ret = dynamic ? only(@variables $name[1:n](t)) :
                                only(@parameters $name[1:n])
            else
                throw("'n' must be greater than or equal to 1")
            end
        else
            V0_val = ustrip(Float64, mV, V0) #FIXME: assumes V0 <: Voltage
            if n == one(n)
                ret = dynamic ? only(@variables $name(t) = V0_val) :
                                only(@parameters $name = V0_val)
            elseif n > one(n)
                ret = dynamic ? only(@variables $name[1:n](t) = V0_val) :
                                only(@parameters $name[1:n] = V0_val)
            else
                throw("'n' must be greater than or equal to 1")
            end
        end
    
        ret = set_symarray_metadata(ret, PrimitiveSource, source)
        ret = set_symarray_metadata(ret, MembranePotential, true)
        return ret
    end
end

ExtrinsicPotential(;V0 = nothing, n = 1, dynamic = true,
                   source::PrimitiveSource = Extrinsic, name::Symbol = :Vₓ) = 
MembranePotential(V0; dynamic = dynamic, source = source, n = n, name = name)

isvoltage(x) = hasmetadata(value(x), MembranePotential)
isintrinsic(x) = hasmetadata(value(x), PrimitiveSource) ?
                 getmetadata(value(x), PrimitiveSource) == Intrinsic : false
isextrinsic(x) = hasmetadata(value(x), PrimitiveSource) ?
                 getmetadata(value(x), PrimitiveSource) == Extrinsic : false


# Custom Unitful.jl quantities
@derived_dimension SpecificConductance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^3 # conductance per unit area
@derived_dimension SpecificCapacitance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^4 # capacitance per unit area
@derived_dimension ConductancePerFarad 𝐓^-1 # S/F cancels out to 1/s; perhaps find a better abstract type?

#const TimeF64 = Quantity{Float64, 𝐓, U} where U

# Ion species
# TODO: Come up with a better solution for this. It should be more easily user extendable
@enum IonSpecies::UInt128 begin
    NonIonic    = 1 << 0
    Sodium      = 1 << 1
    Potassium   = 1 << 2
    Chloride    = 1 << 3
    Calcium     = 1 << 4
    Cation      = 1 << 5
    Anion      = 1 << 6
    Glutamatergic = 1 << 7
    Cholinergic = 1 << 8
    AMPA        = 1 << 9
    NMDA        = 1 << 10
end

const Ca = Calcium
const Na = Sodium
const K = Potassium
const Cl = Chloride
const Mixed = const Leak = NonIonic

const PERIODIC_SYMBOL = IdDict(Na => :Na, K  => :K, Cl => :Cl, Ca => :Ca, Leak => :l)

struct IonConcentration
    ion::IonSpecies
    loc::PrimitiveLocation
end

const Concentration = IonConcentration

function IonConcentration(
    ion::IonSpecies,
    val = nothing;
    location::PrimitiveLocation = Inside,
    dynamic::Bool = false,
    name::Symbol = PERIODIC_SYMBOL[ion]
)

    sym = Symbol(name,(location == Inside ? "ᵢ" : "ₒ"))
    var = dynamic ? only(@variables $sym(t)) : only(@parameters $sym) 
    var = setmetadata(var,  IonConcentration, IonConcentration(ion, location))
    if !isnothing(val)
        if val isa Molarity
            var = setmetadata(var, ConductorUnits, unit(val))
            # FIXME: use proper unit checking
            raw_val = ustrip(Float64, µM, val)
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
    dynamic::Bool = true,
    name::Symbol = Symbol("I", PERIODIC_SYMBOL[ion])
)
    var = dynamic ? only(@variables $name(t)) : only(@parameters $name)
    var = setmetadata(var, IonCurrent, IonCurrent(ion, aggregate))
    if !isnothing(val)
        if val isa Current
            var = setmetadata(var, ConductorUnits, unit(val))
            # FIXME: use proper unit checking
            raw_val = ustrip(Float64, µA, val)
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
iscurrent(x::IonCurrent) = true
getcurrent(x) = iscurrent(x) ? getmetadata(value(x), IonCurrent) : nothing
getcurrent(x::IonCurrent) = x
getion(x::IonCurrent) = getfield(getcurrent(x), :ion)
isaggregate(x::Num) = iscurrent(x) ? getfield(getcurrent(x), :agg) : false

struct EquilibriumPotential
    ion::IonSpecies
end

const Equilibrium = EquilibriumPotential

function EquilibriumPotential(ion::IonSpecies, val; dynamic = false, name::Symbol = PERIODIC_SYMBOL[ion])
    sym = Symbol("E", name)
    var = dynamic ? only(@variables $sym(t)) : only(@parameters $sym) 
    var = setmetadata(var, EquilibriumPotential, EquilibriumPotential(ion))
    if !isnothing(val)
        if val isa Voltage
            var = setmetadata(var, ConductorUnits, unit(val))
            # FIXME: do proper unit checking
            raw_val = ustrip(Float64, mV, val)
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
getion(x::EquilibriumPotential) = getfield(x, :ion)

function getion(x)
    iscurrent(x) && return getion(getcurrent(x))
    isreversal(x) && return getion(getreversal(x))
    return nothing
end

# Alternate constructor; this needs to be better...
function Equilibria(equil::Vector)
    out = Num[]
    for x in equil
        x.first isa IonSpecies || throw("Equilibrium potential must be associated with an ion type.")
        if x.second isa Tuple
            tup = x.second
            tup[2] isa Symbol || throw("Second tuple argument for $(x.first) must be a symbol.")
            push!(out, Equilibrium(x.first, tup[1], dynamic = tup[1] isa Voltage ? false : true, name = tup[2]))
        else
            push!(out, Equilibrium(x.first, x.second, dynamic = x.second isa Voltage ? false : true))
        end
    end
    return out
end
