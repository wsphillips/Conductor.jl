
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

conductance(x::Synapse{ConductanceSystem}) = x.system
reversal(x::Synapse{ConductanceSystem}) = x.metadata[:reversal]

abstract type AbstractJunction end

struct Junction{T} <: AbstractJunction
    system::T
    metadata::Dict{Symbol,Any}
end

function Junction(x::ConductanceSystem, rev::Num)
    metadata = Dict([:reversal => rev]) 
    return Junction{typeof(x)}(x, metadata)
end

conductance(x::Junction{ConductanceSystem}) = x.system
reversal(x::Junction{ConductanceSystem}) = x.metadata[:reversal]

struct Arborization
    parent::Union{Nothing, Junction}
    children::Vector{Junction}
end

Arborization() = Arborization(nothing, Junction[])

# LOCAL / 1st degree connected nodes

function compartments(x::Arborization)
    out = []
    compartments!(out, x)
    return out
end

function compartments!(out::Vector, x::Arborization)
    isnothing(x.parent) || push!(out, x.parent.system)
    append!(out, getproperty.(x.children, :system))
    return
end

function conductances(x::Arborization)
    out = []
    conductances!(out, x)
    return out
end

function conductances!(out::Vector, x::Arborization)
    isnothing(x.parent) || push!(out, x.parent.conductance)
    append!(out, getproperty.(x.children, :conductance))
    return
end

function reversals(x::Arborization)
    out = []
    reversals!(out, x)
    return out
end

function reversals!(out::Vector, x::Arborization)
    comps = compartments(x)
    for comp in comps
        v = get_voltage(comp)
        push!(out, getvar(getname(v)))
    end
    return
end

include("dynamics.jl")

"""
$(TYPEDEF)

A neuronal compartment.

$(TYPEDFIELDS)
"""
struct CompartmentSystem{T <: AbstractDynamics} <: AbstractCompartmentSystem
    "Voltage potential."
    voltage::Num
    dynamics::T
    "Membrane capacitance."
    capacitance::Union{Nothing,Num}
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
    synapses::Vector
    arbor::Arborization
    stimuli::Vector
    "Morphological geometry of the compartment."
    geometry::Geometry
    function CompartmentSystem(voltage, dynamics::T, capacitance, eqs, iv, states, ps,
                               observed, name, systems, defaults, extensions, synapses,
                               arbor, stimuli, geometry;
                               checks = false) where {T <: AbstractDynamics}
        if checks
            foreach(x -> isreversal(x) || throw("Invalid Equilibrium Potential"),
                    channel_reversals)
            foreach(x -> iscurrent(x.lhs) || throw("Invalid current stimulus"), stimuli)
        end
        new{T}(voltage, dynamics, capacitance, eqs, iv, states, ps, observed, name, systems,
               defaults, extensions, synapses, arbor, stimuli, geometry)
    end
end

const Compartment = CompartmentSystem

# can be a `CurrentSystem` constructor
function process_stimulus!(currents, gen, stimulus::PulseTrain)
    (; eqs, dvs, ps, defs) = gen

    I = only(get_variables(stimulus.lhs))
    _vars = filter!(x -> !isequal(x, t), get_variables(stimulus.rhs))
    foreach(x -> push!(isparameter(x) ? ps : dvs, x), _vars)
    hasdefault(I) || push!(defs, I => stimulus.rhs)
    if isparameter(I)
        push!(ps, I)
        push!(currents, -I)
    else
        validate(stimulus) && push!(eqs, stimulus)
        push!(dvs, I)
        push!(currents, -I)
    end
    return
end

# current bias stimulus (can be a `CurrentSystem` constructor)
function process_stimulus!(currents, gen, sym, stimulus::Bias{T}) where {T <: Current}
    ps = gen.ps
    push!(ps, sym)
    push!(currents, -sym)
    return
end

function process_stimulus!(currents, gen, sym, stimulus::PulseTrain{T}) where {T <: Current}
    (; dvs, eqs) = gen
    push!(dvs, sym)
    push!(eqs, sym ~ current_pulses(t, stimulus))
    push!(currents, -sym)
end

function process_stimuli!(currents, gen, stimuli)
    for sym in stimuli
        process_stimulus!(currents, gen, sym, get_stimulus(sym))
    end
end

function resolve_channel_inputs!(gen, dynamics::HodgkinHuxley, synapses)
    all_inputs = Set{Num}()
    # Resolve in/out: "connect" / auto forward cell states to channels
    for chan in union(dynamics.channels, conductance.(synapses))
        for inp in get_inputs(chan)
            isextrinsic(inp) && continue
            push!(all_inputs, inp)
            subinp = getproperty(chan, tosymbol(inp, escape = false))
            push!(gen.eqs, inp ~ subinp)
        end
    end
    return all_inputs
end

function get_required_states!(gen, dynamics::HodgkinHuxley, synapses)
    required_states = resolve_channel_inputs!(gen, dynamics, synapses)
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
    (; eqs, dvs, ps) = gen
    component_currents = filter(x -> iscurrent(x) && !isaggregate(x), union(dvs, ps))
    # Resolve unavailable states
    for s in required_states
        if iscurrent(s) && isaggregate(s)
            eq = s ~ sum(filter(x -> getion(x) == getion(s), component_currents))
            validate(eq) && push!(eqs, eq)
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
function CompartmentSystem(Vₘ::Num, dynamics::T; defaults = Dict(),
                           capacitance = 1µF/cm^2,
                           geometry::Geometry = Point(),
                           stimuli::Vector = [],
                           extensions::Vector{ODESystem} = ODESystem[],
                           name::Symbol = Base.gensym("Compartment"),
                          ) where {T <: AbstractDynamics}

    @parameters cₘ=ustrip(Float64, µF / cm^2, capacitance) [unit = µF / cm^2]
    arbor = Arborization()
    synapses = []
    return CompartmentSystem(Vₘ, dynamics, synapses, arbor, cₘ, geometry, stimuli, defaults,
                             extensions, name)
end

struct LIF <: AbstractDynamics
    V::Num
    τ_mem::Num
    τ_syn::Num
    V_th::Num
    R::Num
    inputs::Vector{Num}
    stimuli::Vector{Num}
end

function LIF(voltage = 0.0; tau_membrane = 10.0, tau_synaptic = 10.0, threshold = 1.0,
             resistance = 1.0, stimulus = 1.0)
    @variables V(t) = voltage
    @parameters V_th = threshold
    @parameters τ_mem=tau_membrane τ_syn=tau_synaptic R=resistance Iₑ=stimulus
    inputs = Num[]
    return LIF(V, τ_mem, τ_syn, V_th, R, inputs, [Iₑ])
end

function CompartmentSystem(dynamics::LIF, defaults, extensions, name)
    (; V, τ_mem, τ_syn, V_th, R, stimuli) = dynamics
    Iₑ = stimuli[1]
    @variables I(t)=0.0 S(t)=false
    @parameters V_rest = MTK.getdefault(V)
    gen = GeneratedCollections(dvs = Set((V, I)),
                               ps = Set((τ_mem, τ_syn, V_th, R, V_rest, Iₑ)),
                               eqs = [
                                   D(V) ~ (-(V - V_rest) / τ_mem) + (R * I + R * Iₑ) / τ_mem,
                                   D(I) ~ -I / τ_syn])

    (; eqs, dvs, ps, observed, systems, defs) = gen
    push!(observed, S ~ V >= V_th)
    merge!(defs, defaults)
    return CompartmentSystem(dynamics, eqs, t, collect(dvs), collect(ps), observed, name,
                             systems, defs, extensions)
end

function Base.convert(::Type{ODESystem}, compartment::CompartmentSystem{LIF};
                      with_cb = false)
    dvs = get_states(compartment)
    ps = get_ps(compartment)
    eqs = get_eqs(compartment)
    defs = get_defaults(compartment)
    obs = get_observed(compartment)
    syss = convert.(ODESystem, get_systems(compartment))
    V = @nonamespace compartment.V
    S = @nonamespace compartment.S
    V_th = @nonamespace compartment.V_th
    V_rest = @nonamespace compartment.V_rest

    cb = with_cb ? MTK.SymbolicDiscreteCallback(V >= V_th, [V ~ V_rest]) :
         MTK.SymbolicDiscreteCallback[]

    return ODESystem(eqs, t, dvs, ps; systems = syss, defaults = defs,
                     name = nameof(compartment), discrete_events = cb,
                     observed = obs, checks = CheckComponents)
end

function CompartmentSystem(Vₘ, dynamics::HodgkinHuxley, synapses, arbor, cₘ,
                           geometry::Geometry, stimuli::Vector, defaults, extensions, name)

    membrane_channels = dynamics.channels
    synaptic_channels = conductance.(synapses)
    axial_conductances = conductances(arbor)
    @parameters aₘ=area(geometry) [unit = cm^2]
    gen = GeneratedCollections(dvs = Set(Vₘ), ps = Set((aₘ, cₘ)),
                               systems = union(membrane_channels, synaptic_channels,
                                               axial_conductances))

    (; eqs, dvs, ps, observed, systems, defs) = gen

    # Compile dynamics into system equations, states, parameters, etc
    process_reversals!(gen, dynamics, synapses, arbor)
    currents = generate_currents!(gen, dynamics, synapses, arbor, Vₘ, aₘ)
    process_stimuli!(currents, gen, stimuli)

    # Voltage equation
    eq = D(Vₘ) ~ -sum(currents) / (cₘ * aₘ)
    validate(eq) && push!(eqs, eq)

    # Apply extensions
    foreach(Base.Fix1(extend!,gen), extensions)

    # Resolve undefined states
    required_states = get_required_states!(gen, dynamics, synapses)
    resolve_states!(gen, required_states)

    merge!(defs, defaults)

    return CompartmentSystem(Vₘ, dynamics, cₘ, eqs, t, collect(dvs), collect(ps), observed,
                             name, systems, defs, extensions, synapses, arbor, stimuli,
                             geometry)
end

function SciMLBase.remake(sys::CompartmentSystem;
                          Vₘ = get_voltage(sys),
                          dynamics = get_dynamics(sys),
                          synapses = get_synapses(sys),
                          arbor = get_arbor(sys),
                          capacitance = get_capacitance(sys),
                          geometry = get_geometry(sys),
                          stimuli = get_stimuli(sys),
                          extensions = get_extensions(sys),
                          name = nameof(sys),
                          defaults = get_defaults(sys))
    CompartmentSystem(Vₘ, dynamics, synapses, arbor, capacitance, geometry, stimuli, defaults,
                      extensions, name)
end

get_dynamics(x::AbstractCompartmentSystem) = getfield(x, :form)
get_extensions(x::AbstractCompartmentSystem) = getfield(x, :extensions)

get_geometry(x::CompartmentSystem) = getfield(x, :geometry)
area(x::CompartmentSystem) = only(@parameters aₘ = area(get_geometry(x)))
get_capacitance(x::CompartmentSystem) = getfield(x, :capacitance)

get_output(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).voltage
get_channels(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).channels
get_synapses(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).synaptic_channels
get_stimuli(x::CompartmentSystem{HodgkinHuxley}) = get_dynamics(x).stimuli

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
Base.iterate(comp::CompartmentSystem, state = 1) = state > 1 ? nothing : (comp, state + 1)
function Base.iterate(rev_comp::Iterators.Reverse{CompartmentSystem}, state = 1)
    state < 1 ? nothing : (rev_comp, state - 1)
end
Base.getindex(comp::CompartmentSystem, i::Int) = i == 1 ? comp : throw(BoundsError(P, i))
Base.firstindex(comp::CompartmentSystem) = 1
Base.lastindex(comp::CompartmentSystem) = 1
Base.getindex(comp::CompartmentSystem, i::Number) = comp[convert(Int, i)]
Base.getindex(comp::CompartmentSystem, I) = [comp[i] for i in I]

function Base.convert(::Type{ODESystem}, compartment::CompartmentSystem)
    dvs = get_states(compartment)
    ps = get_ps(compartment)
    eqs = get_eqs(compartment)
    defs = get_defaults(compartment)
    syss = convert.(ODESystem, get_systems(compartment))

    return ODESystem(eqs, t, dvs, ps; systems = syss, defaults = defs,
                     name = nameof(compartment), checks = CheckComponents)
end
