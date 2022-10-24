
struct CurrentSystem <: AbstractCurrentSystem
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
    inputs::Vector{Num}
    ion::IonSpecies
    aggregate::Bool
    function CurrentSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                           output, inputs, ion, aggregate; checks = false)
        if checks
            # placeholder
        end
        new(eqs, iv, states, ps, observed, name, systems, defaults, output, inputs, ion,
            aggregate)
    end
end

ion(x::CurrentSystem) = getfield(x, :ion)
get_extensions(x::AbstractCurrentSystem) = getfield(x, :extensions)
get_inputs(x::CurrentSystem) = getfield(x, :inputs)
get_output(x::CurrentSystem) = getfield(x, :output)

# Main method for conductance based construction
function CurrentSystem(Vₘ::Num, cond::ConductanceSystem, Erev::Num;
                       aₘ = 1, extensions::Vector{ODESystem} = ODESystem[],
                       defaults = Dict(), name::Symbol = nameof(cond))

    # Extend the conductance system to a current system
    (; eqs, dvs, ps, systems, observed, defs) = copy_collections(cond)
    I = IonCurrent(cond)
    g = get_output(cond)
    eq = I ~ g * (Vₘ - Erev) * (1 * get_unit(g) isa SpecificConductance ? aₘ : 1)
    validate(eq) && push!(eqs, eq)
    push!(dvs, I)
    merge!(defs, defaults)
    return CurrentSystem(eqs, t, dvs, ps, observed, name, systems, defs, I, inputs(cond),
                         permeability(cond), isaggregate(cond); checks = false)
end

function CurrentSystem(Vₘ::Num, syn::AbstractSynapse; kwargs...)
    synaptic_cond, synaptic_rev = conductance(syn), reversal(syn)
    return CurrentSystem(Vₘ, synaptic_cond, synaptic_rev; kwargs...)
end

function CurrentSystem(Vₘ::Num, jxn::AbstractJunction; kwargs...)
    axial_cond, branch_Vm = conductance(jxn), reversal(jxn) 
    return CurrentSystem(Vₘ, axial_cond, branch_Vm; kwargs...)
end

# Specializations for stimuli
function CurrentSystem(Vₘ::Num, stimulus::Bias{T};
                       name::Symbol = Base.gensym("electrode")) where {T <: Current}
    (; eqs, dvs, ps, systems, observed, defs) = GeneratedCollections()
    Iₑ = IonCurrent(stimulus) # creates symbolic with stored default/metadata
    push!(ps, Iₑ)
    return CurrentSystem(eqs, t, dvs, ps, observed, stimulus.name, [], defs, Iₑ, [],
                         NonIonic, false; checks = false)
end

# stimulus (pulse train)
function CurrentSystem(Vₘ::Num, stimulus::PulseTrain{T}) where {T <: Current}
    (; eqs, dvs, ps, systems, observed, defs) = GeneratedCollections()
    Iₑ = IonCurrent(stimulus) # creates symbolic with stored default/metadata
    push!(dvs, Iₑ)
    push!(eqs, Iₑ ~ current_pulses(t, stimulus))
    return CurrentSystem(eqs, t, dvs, ps, observed, stimulus.name, [], defs, Iₑ, [],
                         NonIonic, false; checks = false)
end
