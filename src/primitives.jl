# Custom Unitful.jl quantities
@derived_dimension SpecificConductance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^3
@derived_dimension SpecificCapacitance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^4
@derived_dimension ConductancePerFarad 𝐓^-1 # S/F cancels out to 1/s; perhaps find a better abstract type?

@doc "Conductance per unit area." SpecificConductance 
@doc "Capacitance per unit area." SpecificCapacitance

"""
Faraday's Constant

The electric charge of one mole of electrons.

Unicode ℱ can be typed by writing \\scrF then pressing tab in the Julia REPL, and in many editors.

# Examples
```julia-repl
julia> Conductor.ℱ
96485.33212331001 C mol^-1
```
"""
const ℱ = Unitful.q*Unitful.Na # Faraday's constant

"""
The independent variable for time, ``t``.
"""
const t = let name = :t; only(@variables $name [unit=ms]) end

"""
Differential with respect to time, ``t``.
"""
const D = Differential(t)

# TODO: Make IonSpecies user extendable instead of a fixed set of enums
@enum IonSpecies::UInt128 begin
    NonIonic      = 1 << 0
    Sodium        = 1 << 1
    Potassium     = 1 << 2
    Chloride      = 1 << 3
    Calcium       = 1 << 4
    Cation        = 1 << 5
    Anion         = 1 << 6
    Glutamatergic = 1 << 7
    Cholinergic   = 1 << 8
    AMPA          = 1 << 9
    NMDA          = 1 << 10
end

ion_doc = """
Ion species to annotate `ConductanceSystem`, `IonCurrent`, `IonConcentration`, etc. May be
one of:

"""

for x in instances(IonSpecies)
    global ion_doc *= "- "*string(x)*'\n'
end

@doc ion_doc IonSpecies

const Ca = Calcium
const Na = Sodium
const K = Potassium
const Cl = Chloride
const Mixed = const Leak = NonIonic

const PERIODIC_SYMBOL = IdDict(Na => :Na, K  => :K, Cl => :Cl, Ca => :Ca, Leak => :l)

# Properties
@enum PrimitiveSource Intrinsic Extrinsic
@enum PrimitiveLocation Outside Inside

"""
    MembranePotential(V0 = -60mV; <keyword arguments>)

The voltage in an arbitrary compartment.

If `V0 == nothing`, the default value of the resulting variable will be left unassigned.

# Arguments
- `dynamic::Bool = true`: when false, the voltage will be a static parameter.
- `source::PrimitiveSource = Intrinsic`: the expected origin of a voltage state. `Intrinsic`
  sources are states from the parent compartment. `Extrinsic` sources come from other
  compartments (e.g. presynaptic compartments).
- `n::Integer = 1`: when `n > 1`, the voltage will be a symbolic array of length `n`.
- `name::Symbol = :Vₘ`: the symbol to use for the symbolic variable
"""
struct MembranePotential
    function MembranePotential(V0::Union{Nothing,Real,Voltage} = -60mV; dynamic = true,
                               source::PrimitiveSource = Intrinsic, n::Integer = 1,
                               name::Symbol = :Vₘ)

        V0_val = V0 isa Voltage ? ustrip(Float64, mV, V0) : V0
        if n == one(n) 
            if isnothing(V0)
                ret = only(dynamic ? @variables($name(t),[unit=mV]) : @parameters($name,[unit=mV]))
            else
                ret = only(dynamic ? @variables($name(t)=V0_val,[unit=mV]) :
                                     @parameters($name=V0_val,[unit=mV]))
            end
            ret = setmetadata(ret, PrimitiveSource, source)
            ret = setmetadata(ret, MembranePotential, true)
        elseif n > one(n)
            if isnothing(V0)
                ret = only(dynamic ? @variables($name(t)[1:n],[unit=mV]) : @parameters($name[1:n],[unit=mV]))
            else
                ret = only(dynamic ? @variables($name(t)[1:n] = fill(V0_val, n),[unit=mV]) :
                                     @parameters($name[1:n] = fill(V0_val, n),[unit=mV]))
            end
            ret = set_symarray_metadata(ret, PrimitiveSource, source)
            ret = set_symarray_metadata(ret, MembranePotential, true)
        else
            throw("'n' must be greater than or equal to 1")
        end
        return ret
    end
end

# Internal API: Trait queries for MembranePotential
isvoltage(x)   = hasmetadata(value(x), MembranePotential)
isintrinsic(x) = hasmetadata(value(x), PrimitiveSource) ?
                 getmetadata(value(x), PrimitiveSource) == Intrinsic : false
isextrinsic(x) = hasmetadata(value(x), PrimitiveSource) ?
                 getmetadata(value(x), PrimitiveSource) == Extrinsic : false



function find_voltage(vars, source=isintrinsic)
    idx = findfirst(x -> isvoltage(x) && source(x), vars)
    var = unwrap(vars[idx])
    if operation(var) == getindex && symtype(first(arguments(var))) <: AbstractArray
        return wrap(first(arguments(var)))
    else
        return wrap(var)
    end
end

"""
    ExtrinsicPotential(; n = 1, name::Symbol = :Vₓ)

A voltage derived from an external source (i.e. not the parent compartment).

Equivalent to: `MembranePotential(nothing; dynamic=true, source=Extrinsic, n=n, name=name)`

# Arguments
- `n::Integer = 1`: when `n > 1`, the voltage will be a symbolic array of length `n`.
- `name::Symbol = :Vₓ`: the symbol to use for the symbolic variable
"""
function ExtrinsicPotential(; n = 1, name::Symbol = :Vₓ)
    return MembranePotential(nothing; dynamic = true, source = Extrinsic, n = n, name = name)
end

struct IonConcentration
    ion::IonSpecies
    loc::PrimitiveLocation
end

const Concentration = IonConcentration

"""
    IonConcentration(ion::IonSpecies, val = nothing; <keyword arguments>)

An intra/extracellular concentration of ions.

# Arguments
- `location::PrimitiveLocation = Inside`: location (`Inside` or `Outside`) w.r.t. the parent
  compartment (intracellular or extracellular).
- `dynamic::Bool = false`: when false, the concentration will be a static parameter.
- `name::Symbol = Conductor.PERIODIC_SYMBOL[ion]`: the symbol to use for the symbolic
  variable. By default, a lookup table is used to find the ion's symbol on the periodic
  table of elements.
"""
function IonConcentration(ion::IonSpecies, conc::Union{Nothing,Real,Molarity} = nothing;
                          location::PrimitiveLocation = Inside, dynamic::Bool = false,
                          name::Symbol = PERIODIC_SYMBOL[ion])

    sym = Symbol(name,(location == Inside ? "ᵢ" : "ₒ"))
    conc_val = conc isa Molarity ? ustrip(µM, conc) : conc 
    if isnothing(conc_val)
        ret = dynamic ? only(@variables $sym(t) [unit=µM]) : only(@parameters $sym [unit=µM]) 
    else 
        ret = dynamic ? only(@variables($sym(t)=conc_val,[unit=µM])) :
                        only(@parameters($sym=conc_val,[unit=µM])) 
    end
    ret = setmetadata(ret,  IonConcentration, IonConcentration(ion, location))
    return ret
end

# Internal API: Concentration trait queries
isconc(x) = hasmetadata(value(x), IonConcentration)
getconc(x) = isconc(x) ? getmetadata(value(x), IonConcentration) : nothing

struct IonCurrent
    ion::IonSpecies
    agg::Bool
end

"""
    IonCurrent(ion::IonSpecies, val = nothing; <keyword arguments>)

An ionic membrane current.

# Arguments
- `aggregate::Bool = false`: aggregate currents are the sum of all conductances (with 
  matched ion species) flowing into the parent compartment. For example, an aggregate 
  `IonCurrent` for Calcium will be the sum of all other Calcium-permeable currents.
- `dynamic::Bool = true`: when `dynamic == false` the `IonCurrent` will be a static 
  parameter value.
- `name::Symbol = Symbol("I", Conductor.PERIODIC_SYMBOL[ion])`: the symbol to use for the 
  symbolic variable. By default, a lookup table is used to find the ion's symbol on the 
  periodic table of elements.
"""
function IonCurrent(ion::IonSpecies, curr::Union{Nothing,Real,Current} = nothing;
                    aggregate::Bool = false, dynamic::Bool = true,
                    name::Symbol = Symbol("I", PERIODIC_SYMBOL[ion]))

    curr_val = curr isa Current ? ustrip(µA, curr) : curr 

    if isnothing(curr_val)
        ret = dynamic ? only(@variables $name(t) [unit=µA]) : only(@parameters $name [unit=µA]) 
    else
        curr_unit = curr isa Current ? µA : unit(curr)
        ret = dynamic ? only(@variables($name(t) = curr_val,[unit=curr_unit])) :
                        only(@parameters($name = curr_val,[unit=curr_unit])) 
    end
    return setmetadata(ret, IonCurrent, IonCurrent(ion, aggregate))
end

IonCurrent(cond::AbstractConductanceSystem) = 
IonCurrent(permeability(cond); name = Symbol("I", nameof(cond)))

# Internal API: Current trait queries
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

"""
    EquilibriumPotential(ion::IonSpecies, val; <keyword arguments>)

An equilibrium (a.k.a. reversal) potential.

# Arguments
- `dynamic::Bool = false`: a dynamic `EquilbriumPotential` is assumed to vary with time
  (e.g. derived from the Nernst equation).
- `name::Symbol = Symbol("I", Conductor.PERIODIC_SYMBOL[ion])`: the symbol to use for the 
  symbolic variable. By default, a lookup table is used to find the ion's symbol on the 
  periodic table of elements.
"""
function EquilibriumPotential(ion::IonSpecies, eqv::Union{Nothing,Real,Voltage}; dynamic::Bool = false,
                              name::Symbol = PERIODIC_SYMBOL[ion])
    sym = Symbol("E", name)
    eqv_val = eqv isa Voltage ? ustrip(mV, eqv) : eqv 
    if isnothing(eqv_val)
        ret = dynamic ? only(@variables $sym(t) [unit=mV]) : only(@parameters $sym [unit=mV]) 
    else 
        ret = dynamic ? only(@variables($sym(t) = eqv_val,[unit=mV])) :
                        only(@parameters($sym = eqv_val,[unit=mV]))
    end
    return setmetadata(ret, EquilibriumPotential, EquilibriumPotential(ion))
end

# Internal API: EquilibriumPotential trait queries
isreversal(x) = hasmetadata(value(x), EquilibriumPotential)
getreversal(x) = isreversal(x) ? getmetadata(value(x), EquilibriumPotential) : nothing
getion(x::EquilibriumPotential) = getfield(x, :ion)

function getion(x)
    iscurrent(x) && return getion(getcurrent(x))
    isreversal(x) && return getion(getreversal(x))
    return nothing
end

function find_reversal(cond::AbstractConductanceSystem, reversals)
    ion = permeability(cond)
    return only(filter(x -> isequal(getion(x), ion), reversals))
end

#FIXME: this is a kludge
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
