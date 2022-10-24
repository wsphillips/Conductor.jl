
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

# Membrane and reversal potential source from compartment system
function CurrentSystem(Vₘ::Num, cond::ConductanceSystem, Erev::Num;
                       aₘ = 1, extensions::Vector{ODESystem} = ODESystem[],
                       defaults = Dict(), name::Symbol = nameof(cond))

    # Extend the conductance system to a current system
    gen = GeneratedCollections(eqs = get_eqs(cond),
                               systems = get_systems(cond),
                               observed = get_observed(cond),
                               dvs = get_dvs(cond),
                               ps = get_ps(cond),
                               defs = get_defaults(cond))
    (; eqs, dvs, ps, systems, observed, defs) = gen
    I = IonCurrent(cond)
    g = renamespace(cond, get_output(cond))
    eq = I ~ g * (Vₘ - Erev) * (1 * get_unit(g) isa SpecificConductance ? aₘ : 1)
    validate(eq) && push!(eqs, eq)
    push!(dvs, I)
    merge!(defs, defaults)
    return CurrentSystem(eqs, t, dvs, ps, observed, name, systems, defs, I, inputs(cond),
                         conductance(cond), isaggregate(cond); checks = false)
end

# synapse
function CurrentSystem()
    
    return CurrentSystem(eqs, t, dvs, ps, observed, name, systems, defs, I, inputs(cond),
                         conductance(cond), isaggregate(cond); checks = false)
end

# axial
function CurrentSystem()
    
    return CurrentSystem(eqs, t, dvs, ps, observed, name, systems, defs, I, inputs(cond),
                         conductance(cond), isaggregate(cond); checks = false)
end

# stimulus (bias)
function CurrentSystem(sym, stimulus::Bias{T} = get_stimulus(sym)) where {T <: Current}
    ps = Num[]
    outputs = Num[]
    push!(ps, sym)
    push!(currents, -sym)

    return CurrentSystem(eqs, t, dvs, ps, observed, name, systems, defs, I, inputs(cond),
                         conductance(cond), isaggregate(cond); checks = false)
end

# stimulus (pulse train)
function CurrentSystem(stimulus::PulseTrain{T}) where {T <: Current}
    
    return CurrentSystem(eqs, t, dvs, ps, observed, name, systems, defs, I, inputs(cond),
                         conductance(cond), isaggregate(cond); checks = false)
end
