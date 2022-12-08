# Current specifications
abstract type AbstractSynapse end

struct Synapse{T} <: AbstractSynapse
    system::T
    metadata::Dict{Symbol,Any}
end

#TODO: relax these fields to a current ?? 
function Synapse(x::ConductanceSystem, rev::Num)
    metadata = Dict([:reversal => rev])
    return Synapse{typeof(x)}(x, metadata)
end

conductance(x::Synapse) = x.system
reversal(x::Synapse) = x.metadata[:reversal]
conductances(x::Vector{<:Synapse}) = conductance.(x)
reversals(x::Vector{<:Synapse}) = reversal.(x)

abstract type AbstractJunction end

struct Junction{T} <: AbstractJunction
    system::T
    metadata::Dict{Symbol,Any}
end

function Junction(x::ConductanceSystem, rev::Num)
    metadata = Dict([:reversal => rev]) 
    return Junction{typeof(x)}(x, metadata)
end

conductance(x::Junction) = x.system
reversal(x::Junction) = x.metadata[:reversal]

struct Arborization
    parent::Union{Nothing, Junction}
    children::Vector{Junction}
end

Arborization() = Arborization(nothing, Junction[])

# LOCAL / 1st degree connected nodes
function conductances(x::Arborization)
    out = []
    conductances!(out, x)
    return out
end

function conductances!(out::Vector, x::Arborization)
    isnothing(x.parent) || push!(out, conductance(x.parent))
    append!(out, conductance.(x.children))
    return
end

function reversals(x::Arborization)
    out = []
    reversals!(out, x)
    return out
end

function reversals!(out::Vector, x::Arborization)
    isnothing(x.parent) || push!(out, reversal(x.parent))
    append!(out, reversal.(x.children))
    return
end

# CurrentSystem + methods
struct CurrentSystem{T<:ConductanceModel} <: AbstractCurrentSystem
    eqs::Vector{Equation}
    "Independent variable. Defaults to time, ``t``."
    iv::Num
    states::Vector{Num}
    ps::Vector{Num}
    observed::Vector{Equation}
    name::Symbol
    systems::Vector{AbstractTimeDependentSystem}
    defaults::Dict
    # Conductor fields
    "Conductance, ``g``, of the system."
    output::Num
    inputs::Set{Num}
    ion::IonSpecies
    model::T
    function CurrentSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                           output, inputs, ion, model; checks = false)
        if checks
            # placeholder
        end
        new{typeof(model)}(eqs, iv, states, ps, observed, name, systems, defaults, output,
                           inputs, ion, model)
    end
end

ion(x::CurrentSystem) = getfield(x, :ion)
get_extensions(x::AbstractCurrentSystem) = getfield(x, :extensions)
get_inputs(x::CurrentSystem) = getfield(x, :inputs)
get_output(x::CurrentSystem) = getfield(x, :output)
get_model(x::CurrentSystem{<:ConductanceModel}) = getfield(x, :model)
output(x::CurrentSystem) = renamespace(x, get_output(x))
inputs(x::CurrentSystem) = renamespace.(x, get_inputs(x))

# Main method for conductance based construction
function CurrentSystem(Vₘ::Num, cond::ConductanceSystem{T}, Erev::Num;
                       aₘ = 1, extensions::Vector{ODESystem} = ODESystem[],
                       defaults = Dict(), name::Symbol = nameof(cond)) where {T}
    
    # Extend the conductance system to a current system
    (; eqs, dvs, ps, systems, observed, defs) = copy_collections(cond)
    push!(ps, aₘ)
    push!(isparameter(Erev) ? ps : dvs, Erev)
    inps = Set(get_inputs(cond))
    push!(inps, Vₘ)
    I = IonCurrent(cond)
    g = get_output(cond)
    eq = I ~ g * (Vₘ - Erev) * (1 * get_unit(g) isa SpecificConductance ? aₘ : 1)
    if T <: IntegratedSynapse
        @parameters W # NOTE: we do NOT provide default weight values here; do it post hoc 
        push!(ps, W)
        eq = eq.lhs ~ eq.rhs * W
    end
    validate(eq) && push!(eqs, eq)
    #push!(dvs, Vₘ, I)
    push!(dvs, I)
    merge!(defs, defaults)
    return CurrentSystem(eqs, t, collect(dvs), collect(ps), observed, name, systems, defs,
                         I, inps, permeability(cond), get_model(cond); checks = false)
end

# Specializations for stimuli
function CurrentSystem(Vₘ::Num, stimulus::Bias{T};
                       name::Symbol = Base.gensym("electrode")) where {T <: Current}
    (; eqs, dvs, ps, systems, observed, defs) = GeneratedCollections()
    Iₑ = IonCurrent(stimulus) # creates symbolic with stored default/metadata
    push!(ps, Iₑ)
    return CurrentSystem(eqs, t, collect(dvs), collect(ps), observed, stimulus.name, [],
                         defs, Iₑ, Set{Num}(), NonIonic, stimulus; checks = false)
end

# stimulus (pulse train)
function CurrentSystem(Vₘ::Num, stimulus::PulseTrain{T}) where {T <: Current}
    (; eqs, dvs, ps, systems, observed, defs) = GeneratedCollections()
    Iₑ = IonCurrent(stimulus) # creates symbolic with stored default/metadata
    push!(dvs, Iₑ)
    push!(eqs, Iₑ ~ current_pulses(t, stimulus))
    return CurrentSystem(eqs, t, collect(dvs), collect(ps), observed, stimulus.name, [],
                         defs, Iₑ, Set{Num}(), NonIonic, stimulus; checks = false)
end

function Base.convert(::Type{ODESystem}, currsys::CurrentSystem)
    dvs = states(currsys)
    ps = parameters(currsys)
    eqs = equations(currsys)
    defs = get_defaults(currsys)
    sys = ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(currsys),
                    checks = CheckComponents)
    #return extend(sys, get_extensions(sys))
    return sys
end

function Base.:(==)(sys1::CurrentSystem, sys2::CurrentSystem)
    sys1 === sys2 && return true
    iv1 = get_iv(sys1)
    iv2 = get_iv(sys2)
    isequal(iv1, iv2) &&
        isequal(nameof(sys1), nameof(sys2)) &&
        _eq_unordered(get_eqs(sys1), get_eqs(sys2)) &&
        _eq_unordered(get_states(sys1), get_states(sys2)) &&
        _eq_unordered(get_ps(sys1), get_ps(sys2)) &&
        all(s1 == s2 for (s1, s2) in zip(get_systems(sys1), get_systems(sys2)))
end


