abstract type CompartmentForm end

struct HodgkinHuxley <: CompartmentForm
    "Voltage potential."
    voltage::Num
    "Membrane capacitance."
    capacitance::Num
    "Morphological geometry of the compartment."
    geometry::Geometry
    "Ionic conductances."
    channels::Vector{AbstractConductanceSystem}
    "Equilibrium potentials belonging to ionic membrane conductances."
    channel_reversals::Vector{Num}
    "Synaptic conductances."
    synaptic_channels::Vector{AbstractConductanceSystem}
    "Equilibrium potentials belonging to synaptic conductances."
    synaptic_reversals::Vector{Num}
    "Axial (intercompartmental) conductances."
    axial_conductances::Vector{Tuple{AbstractConductanceSystem,Num}}
    "Experimental stimuli (for example, current injection)."
    stimuli::Vector{Equation}
end

# basic constructor
function HodgkinHuxley(
    Vₘ::Num,
    channels,
    channel_reversals;
    capacitance = 1µF/cm^2,
    geometry::Geometry = Point(),
    stimuli::Vector{Equation} = Equation[])

    @parameters cₘ = ustrip(Float64, mF/cm^2, capacitance)

    synaptic_channels = Vector{AbstractConductanceSystem}()
    synaptic_reversals = Vector{Num}()
    axial_conductances = Vector{Tuple{AbstractConductanceSystem,Num}}()
   
    HodgkinHuxley(Vₘ, cₘ, geometry, channels, channel_reversals, synaptic_channels,
                  synaptic_reversals, axial_conductances, stimuli)
end

"""
$(TYPEDEF)

A neuronal compartment.

$(TYPEDFIELDS)
"""
struct CompartmentSystem{T<:CompartmentForm} <: AbstractCompartmentSystem
    form::T
    eqs::Vector{Equation}
    "Independent variable. Defaults to time, ``t``."
    iv::Num
    states::Vector{Num}
    ps::Vector{Num}
    observed::Vector{Equation}
    name::Symbol
    systems::Vector{AbstractTimeDependentSystem}
    defaults::Dict
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    "Reference to the parent system when the compartment is a subcompartment in a `MultiCompartmentSystem`."
    parent::Ref{AbstractCompartmentSystem}

    function CompartmentSystem(form::T, eqs, iv, states, ps, observed, name, systems,
                               defaults, extensions, parent; checks = false) where {T <: CompartmentForm}
        if checks
        # placeholder
        foreach(x -> isreversal(x) || throw("Invalid Equilibrium Potential"), channel_reversals)
        foreach(x -> iscurrent(x.lhs) || throw("Invalid current stimulus"), stimuli)

        end
        new{T}(form, eqs, iv, states, ps, observed, name, systems, defaults, extensions, parent)
    end
end

const Compartment = CompartmentSystem

struct GeneratedCollections
    eqs::Vector{Equation}
    dvs::Set{Num}
    ps::Set{Num}
    systems::Vector{AbstractTimeDependentSystem}
    observed::Vector{Equation}
    defs::Dict
end

function GeneratedCollections(; eqs = Equation[], systems = AbstractTimeDependentSystem[],
                              observed = Equation[], dvs = Set{Num}(), ps  = Set{Num}(),
                              defs = Dict())
    return GeneratedCollections(eqs, dvs, ps, systems, observed, defs)
end

function process_reversals!(gen, dynamics)
    reversal_equation_vars = Set{Num}()
    @unpack channel_reversals, synaptic_reversals, axial_conductances = dynamics
    @unpack dvs, ps, eqs = gen
    all_reversals = union(channel_reversals, synaptic_reversals, last.(axial_conductances))

    for Erev in all_reversals
        if isparameter(Erev)
            push!(ps, Erev)
        else
            push!(dvs, Erev)
            # FIXME: hacked solution. This assumes the equation for a dynamic reversal == default value
            if MTK.hasdefault(Erev)
                get_variables!(reversal_equation_vars, getdefault(Erev))
                push!(eqs, Erev ~ getdefault(Erev))
            end
        end
    end
    filter!(x -> !isequal(x, t), reversal_equation_vars) # remove iv
    foreach(x -> isparameter(x) && push!(ps, x), reversal_equation_vars)
    return nothing
end

function generate_currents!(gen, dynamics, Vₘ, aₘ)
    @unpack channels, channel_reversals, synaptic_channels, synaptic_reversals,
            axial_conductances = dynamics
    @unpack eqs, dvs = gen
    paired_channels = zip(channels,
                          broadcast(chan -> find_reversal(chan, channel_reversals), channels))
    paired_synapses = zip(synaptic_channels,
                          broadcast(chan -> find_reversal(chan, synaptic_reversals), synaptic_channels))
    paired_conductances = union(axial_conductances, paired_channels, paired_synapses)
    currents = Set()
    for (chan, Erev) in paired_conductances
        I = IonCurrent(chan)
        g = renamespace(chan, get_output(chan))
        push!(eqs, I ~ g*(Vₘ - Erev)*(1*getmetadata(g, ConductorUnits) isa SpecificConductance ? aₘ : 1))
        push!(dvs, I)
        push!(currents, I)
    end
    return currents
end

function process_stimuli!(currents, gen, dynamics)
    @unpack eqs, dvs, ps, defs = gen
    for stimulus in dynamics.stimuli
        I = only(get_variables(stimulus.lhs))
        _vars = filter!(x -> !isequal(x,t), get_variables(stimulus.rhs))
        foreach(x -> push!(isparameter(x) ? ps : dvs, x), _vars)
        hasdefault(I) || push!(defs, I => stimulus.rhs)
        if isparameter(I)
            push!(ps, I)
            push!(currents, -I)
        else
            push!(eqs, stimulus)
            push!(dvs, I)
            push!(currents, -I)
        end
    end
end

function resolve_channel_inputs!(gen, dynamics)
    all_inputs = Set{Num}()
    # Resolve in/out: "connect" / auto forward cell states to channels
    for chan in union(dynamics.channels, dynamics.synaptic_channels)
        for inp in get_inputs(chan)
            isextrinsic(inp) && continue
            push!(all_inputs, inp)
            subinp = getproperty(chan, tosymbol(inp, escape=false))
            push!(gen.eqs, inp ~ subinp)
        end
    end
    return all_inputs
end

function get_required_states!(gen, dynamics)
    required_states = resolve_channel_inputs!(gen, dynamics)
    internally_modified = Set{Num}()
    foreach(x -> get_variables!(internally_modified, x.lhs), gen.eqs)
    union!(required_states, setdiff(gen.dvs, internally_modified))
    return required_states
end

function extend!(gen, extension)
    union!(gen.eqs, equations(extension))
    union!(gen.ps, parameters(extension))
    union!(gen.dvs, states(extension))
    return nothing
end

function resolve_states!(gen, required_states)
    @unpack eqs, dvs, ps = gen
    component_currents = filter(x -> iscurrent(x) && !isaggregate(x), union(dvs,ps))
    # Resolve unavailable states
    for s in required_states
        if iscurrent(s) && isaggregate(s)
            push!(eqs, s ~ sum(filter(x -> getion(x) == getion(s), component_currents)))
            push!(dvs, s)
        end
        # etc for more switch cases...
    end
    return nothing
end

"""
    CompartmentSystem(Vₘ, channels, reversals; <keyword arguments>)

# Arguments
- `capacitance::SpecificCapacitance`: The capacitance of the compartment given in Farads 
  per unit area (e.g. µF/cm^2).
- `geometry::Geometry`: Morphological geometry of the compartment.
- `extensions::Vector{ODESystem}`: Additional systems to extend dynamics. Extensions are
  composed with the parent system during conversion to `ODESystem`.
- `stimuli::Vector{Equation}`: 
- `name::Symbol`: Name of the system.
"""
function CompartmentSystem(dynamics; defaults = Dict(),
                           extensions::Vector{ODESystem} = ODESystem[],
                           name::Symbol = Base.gensym("Compartment")) 
    parent = Ref{AbstractCompartmentSystem}()
    return CompartmentSystem(dynamics, defaults, extensions, name, parent)
end

struct LIF <: CompartmentForm
    V::Num
    τ_mem::Num
    τ_syn::Num
    V_th::Num
    R::Num
    inputs::Vector{Num}
    stimuli::Vector{Num}
end

function LIF(voltage = 0.0; tau_membrane = 10.0, tau_synaptic = 10.0, threshold = 1.0, resistance = 1.0, stimulus = 1.0)
    @variables V(t) = voltage
    @parameters V_th = threshold
    @parameters τ_mem = tau_membrane τ_syn = tau_synaptic R = resistance Iₑ = stimulus
    inputs = Num[]
    return LIF(V, τ_mem, τ_syn, V_th, R, inputs, [Iₑ])
end

function CompartmentSystem(dynamics::LIF, defaults, extensions, name, parent)
    @unpack V, τ_mem, τ_syn, V_th, R, inputs, stimuli = dynamics
    Iₑ = stimuli[1]
    @variables I(t) = 0.0 S(t) = false
    @parameters V_rest = MTK.getdefault(V)
    gen = GeneratedCollections(dvs = Set((V, I, S)),
                               ps = Set((τ_mem, τ_syn, V_th, R, V_rest, Iₑ)),
                               eqs = [D(V) ~ (-(V-V_rest)/τ_mem) + (R*(I + Iₑ))/τ_mem,
                                      S ~ V == V_th,
                                      D(I) ~ -I/τ_syn + sum(inputs)])
    @unpack eqs, dvs, ps, observed, systems, defs = gen
    merge!(defs, defaults)
    return CompartmentSystem(dynamics, eqs, t, collect(dvs), collect(ps), observed, name,
                             systems, defs, extensions, parent)
end

function CompartmentSystem(dynamics::HodgkinHuxley, defaults, extensions, name, parent)

    @unpack channels, synaptic_channels, axial_conductances, stimuli = dynamics
    @parameters aₘ = area(dynamics.geometry)
    Vₘ, cₘ = dynamics.voltage, dynamics.capacitance
    gen = GeneratedCollections(dvs = Set(Vₘ), ps = Set((aₘ,cₘ)),
                               systems = union(channels, synaptic_channels,
                                               first.(axial_conductances)))
    @unpack eqs, dvs, ps, observed, systems, defs = gen
    
    # Compile dynamics into system equations, states, parameters, etc
    process_reversals!(gen, dynamics)
    currents = generate_currents!(gen, dynamics, Vₘ, aₘ)
    process_stimuli!(currents, gen, dynamics)
    # Voltage equation
    push!(eqs, D(Vₘ) ~ -sum(currents)/(cₘ*aₘ))
    # Apply extensions
    extend_gen! = Base.Fix1(extend!, gen)
    foreach(extend_gen!, extensions)
    # Resolve undefined states
    required_states = get_required_states!(gen, dynamics)
    resolve_states!(gen, required_states)

    merge!(defs, defaults)

    return  CompartmentSystem(dynamics, eqs, t, collect(dvs), collect(ps), observed, name,
                              systems, defs, extensions, parent)
end

function SciMLBase.remake(sys::CompartmentSystem;
                          dynamics = get_dynamics(sys),
                          extensions = get_extensions(sys),
                          parent = parent(sys),
                          name = nameof(sys),
                          defaults = get_defaults(sys))

    CompartmentSystem(dynamics, defaults, extensions, name, parent)
end

get_dynamics(x::AbstractCompartmentSystem) = getfield(x, :form)
get_extensions(x::AbstractCompartmentSystem) = getfield(x, :extensions)

get_geometry(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).geometry
area(x::CompartmentSystem{HodgkinHuxley}) = only(@parameters aₘ = area(get_geometry(x)))
capacitance(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).capacitance

get_output(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).voltage
get_channels(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).channels
get_synapses(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).synaptic_channels
get_stimuli(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).stimuli

hasparent(x::CompartmentSystem) = isassigned(getfield(x, :parent))
parent(x::CompartmentSystem) = getfield(x, :parent)

function setparent!(child::CompartmentSystem, parent::AbstractCompartmentSystem)
    ref = getfield(child, :parent)
    ref[] = parent
    return nothing
end

function get_reversals(x::CompartmentSystem{HodgkinHuxley})
    dyn = get_dynamics(x)
    return dyn.channel_reversals, dyn.synaptic_reversals
end

get_synaptic_reversals(x::CompartmentSystem{HodgkinHuxley}) = get_reversals(x)[2]
get_channel_reversals(x::CompartmentSystem{HodgkinHuxley}) = get_reversals(x)[1]

function Base.:(==)(sys1::AbstractCompartmentSystem, sys2::AbstractCompartmentSystem)
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

Base.eltype(::CompartmentSystem) = CompartmentSystem
Base.length(::CompartmentSystem) = 1
Base.iterate(comp::CompartmentSystem, state=1) = state > 1 ? nothing : (comp, state+1)
Base.iterate(rev_comp::Iterators.Reverse{CompartmentSystem}, state=1) = state < 1 ? nothing : (rev_comp, state-1)
Base.getindex(comp::CompartmentSystem, i::Int) = i == 1 ? comp : throw(BoundsError(P, i))
Base.firstindex(comp::CompartmentSystem) = 1
Base.lastindex(comp::CompartmentSystem) = 1
Base.getindex(comp::CompartmentSystem, i::Number) = comp[convert(Int, i)]
Base.getindex(comp::CompartmentSystem, I) = [comp[i] for i in I]

function Base.convert(::Type{ODESystem}, compartment::CompartmentSystem)

    dvs = get_states(compartment)
    ps  = get_ps(compartment)
    eqs = get_eqs(compartment)
    defs = get_defaults(compartment)
    syss = convert.(ODESystem, get_systems(compartment))
    
    return ODESystem(eqs, t, dvs, ps; systems = syss, defaults = defs,
                     name = nameof(compartment))
end

function Base.convert(::Type{ODESystem}, compartment::CompartmentSystem{LIF})

    dvs = get_states(compartment)
    ps  = get_ps(compartment)
    eqs = get_eqs(compartment)
    defs = get_defaults(compartment)
    syss = convert.(ODESystem, get_systems(compartment))
    
    root_eqs = [@nonamespace(compartment.V) ~ @nonamespace(compartment.V_th)]
    affect = [@nonamespace(compartment.V) ~ @nonamespace(compartment.V_rest)]

    return ODESystem(eqs, t, dvs, ps; systems = syss, defaults = defs,
                     name = nameof(compartment), continuous_events = root_eqs => affect)
end

