
abstract type Geometry end

struct Sphere <: Geometry
    radius
end

function Sphere(; radius)
    Sphere(radius)
end

struct Point <: Geometry end

struct Cylinder <: Geometry
    radius
    height
    open_ends::Bool
end

function Cylinder(; radius, height, open_ends = true)
    return Cylinder(radius, height, open_ends)
end

height(x::Geometry) = isdefined(x, :height) ? getfield(x, :height) : nothing
radius(x::Geometry) = getfield(x, :radius)
radius(::Point) = 0.0

# TODO: Remove unit stripping when we have proper unit checking implemented
area(x::Sphere) = ustrip(Float64, cm^2, 4*π*radius(x)^2)
area(x::Cylinder) = ustrip(Float64, cm^2, 2*π*radius(x)*(height(x) + (x.open_ends ? 0µm : radius(x))))
area(::Point) = 1.0

abstract type AbstractCompartmentSystem <: AbstractTimeDependentSystem end

struct CompartmentSystem <: AbstractCompartmentSystem
    iv::Num # usually just t
    ## Intrinsic properties
    voltage::Num # symbol that represents membrane voltage
    capacitance::Num # specific membrane capacitance
    geometry::Geometry # compartment geometry metadata (shape, dimensions, etc)
    ## Dynamics
    chans::Set{AbstractConductanceSystem} # conductance systems
    channel_reversals::Set{Num}
    synapses::Set{AbstractConductanceSystem} # synaptic conductance systems
    synaptic_reversals::Set{Num}
    axial_conductance::Set{Tuple{AbstractConductanceSystem,Num}}
    stimuli::Vector{Equation}
    extensions::Vector{ODESystem}
    defaults::Dict
    name::Symbol
    eqs::Vector{Equation}
    systems::Vector{AbstractTimeDependentSystem}
    observed::Vector{Equation}
    function CompartmentSystem(iv, voltage, capacitance, geometry, chans, channel_reversals,
                               synapses, synaptic_reversals, axial_conductance, stimuli, extensions, defaults,
                               name, eqs, systems, observed; checks = false)
        if checks
        # placeholder
        end
        new(iv, voltage, capacitance, geometry, chans, channel_reversals, synapses,
            synaptic_reversals, axial_conductance, stimuli, extensions, defaults, name, eqs, systems, observed)
    end
end

const Compartment = CompartmentSystem

function CompartmentSystem(
    Vₘ::Num,
    channels,
    reversals;
    capacitance = 1µF/cm^2,
    geometry::Geometry = Point(),
    extensions::Vector{ODESystem} = ODESystem[],
    stimuli::Vector{Equation} = Equation[],
    name::Symbol = Base.gensym("Compartment")
) 
    @parameters cₘ = ustrip(Float64, mF/cm^2, capacitance)
    foreach(x -> isreversal(x) || throw("Invalid Equilibrium Potential"), reversals)
    foreach(x -> iscurrent(x.lhs) || throw("Invalid current stimulus"), stimuli)
    eqs = Equation[]
    systems = AbstractSystem[]
    observed = Equation[]
    return CompartmentSystem(t, Vₘ, cₘ, geometry, Set(channels), Set(reversals), Set(),
                             Set(), Set(), stimuli, extensions, Dict(), name, eqs, systems, observed)
end

# AbstractSystem interface extensions
get_geometry(x::AbstractCompartmentSystem) = getfield(x, :geometry)
area(x::AbstractCompartmentSystem) = only(@parameters aₘ = area(get_geometry(x)))
capacitance(x::AbstractCompartmentSystem) = getfield(x, :capacitance)
get_output(x::AbstractCompartmentSystem) = getfield(x, :voltage)

# TODO: Implement these...
#get_observed(x::CompartmentSystem) = 0
#get_continuous_events(x::CompartmentSystem) = 0

get_extensions(x::AbstractCompartmentSystem) = getfield(x, :extensions)
get_channels(x::AbstractCompartmentSystem) = getfield(x, :chans)
get_synapses(x::AbstractCompartmentSystem) = getfield(x, :synapses)
get_stimuli(x::AbstractCompartmentSystem) = getfield(x, :stimuli)

# TODO: define top-level and recursive getters for inputs
function get_inputs(x::AbstractCompartmentSystem) end 
function inputs(x::AbstractCompartmentSystem) end

function get_reversals(x::AbstractCompartmentSystem)
    return getfield(x, :channel_reversals), getfield(x, :synaptic_reversals) 
end

get_synaptic_reversals(x::AbstractCompartmentSystem) = get_reversals(x)[2]
get_channel_reversals(x::AbstractCompartmentSystem) = get_reversals(x)[1]

function build_toplevel(system)
    dvs = Set{Num}()
    ps = Set{Num}()
    eqs = Equation[]
    defs = Dict()
    build_toplevel!(dvs, ps, eqs, defs, system)
end

function build_toplevel!(dvs, ps, eqs, defs, comp_sys::CompartmentSystem)

    currents = Set{Num}()
    syncurrents = Set{Num}()
    aₘ = area(comp_sys)
    cₘ = capacitance(comp_sys)
    Vₘ = get_output(comp_sys)
    push!(dvs, Vₘ)
    push!(ps, aₘ)
    push!(ps, cₘ)
    chan_revs, syn_revs = get_reversals(comp_sys) 

    # Check for possible dynamic reversals
    rev_eq_vars = []
    for Erev in union(chan_revs, syn_revs)
        if isparameter(Erev)
            push!(ps, Erev)
        else
            push!(dvs, Erev)
            get_variables!(rev_eq_vars, getdefault(Erev))
            push!(eqs, Erev ~ getdefault(Erev))
        end
    end
    
    foreach(x -> isparameter(x) && push!(ps, x), rev_eq_vars)

    # Gather channel current equations
    for chan in get_channels(comp_sys)
        ion = permeability(chan)
        Erev = only(filter(x -> isequal(getion(x), ion), chan_revs))
        I = IonCurrent(ion, name = Symbol("I", nameof(chan)))
        g = renamespace(chan, get_output(chan))
        # TODO: checks for whether we need area scaling
        push!(eqs, I ~ g*aₘ*(Vₘ - Erev))
        push!(dvs, I)
        push!(currents, I)
        merge!(defs, defaults(chan))
    end

    # Gather synaptic current equations
    for synapse in get_synapses(comp_sys)
        ion = permeability(synapse)
        Esyn = only(filter(x -> isequal(getion(x), ion), syn_revs))
        I = IonCurrent(ion, name = Symbol("I", nameof(synapse)))
        s = renamespace(synapse, get_output(synapse))
        push!(eqs, I ~ s*(Vₘ - Esyn))
        push!(dvs, I)
        push!(syncurrents, I)
        merge!(defs, defaults(synapse))
    end
    
    # Gather axial current equations
    for (ax, childvm) in get_axial_conductance(comp_sys)
        I = IonCurrent(Leak, name = Symbol("I", nameof(ax)))
        g = renamespace(ax, get_output(ax))
        # TODO: we should have a metadata check for area-dependent units
        push!(eqs, I ~ g*aₘ*(childvm - Vₘ))
        push!(dvs, I)
        push!(dvs, childvm)
        push!(currents, I)
        merge!(defs, defaults(ax))
    end

    # Gather extension equations
    for extension in get_extensions(comp_sys)
        union!(eqs, equations(extension))
        union!(ps, parameters(extension))
        union!(dvs, states(extension))
        merge!(defs, defaults(extension))
    end
    
    for stimulus in get_stimuli(comp_sys)
        I = only(get_variables(stimulus.lhs))
        hasdefault(I) || push!(defs, I => stimulus.rhs)
        if isparameter(I)
            push!(ps, I)
        else
            push!(eqs, stimulus)
            push!(dvs, I)
            push!(currents, -I)
        end
    end
    # voltage equation
    push!(eqs, D(get_output(comp_sys)) ~ -sum(union(currents,syncurrents))/(cₘ*aₘ))

    return eqs, dvs, ps, defs, currents, syncurrents
end

# collect _top level_ eqs including from extension + currents + reversals + Vₘ
function get_eqs(x::AbstractCompartmentSystem; rebuild = false)
    #if rebuild || isempty(getfield(x, :eqs))
        empty!(getfield(x, :eqs))
        union!(getfield(x, :eqs), build_toplevel(x)[1])
    #end
    return getfield(x, :eqs)
end

function get_states(x::AbstractCompartmentSystem)
    collect(build_toplevel(x)[2])
end

MTK.has_ps(x::CompartmentSystem) = true

function get_ps(x::AbstractCompartmentSystem)
    # collect _top level_ parameters from extension + currents + capacitance + area + reversals
    collect(build_toplevel(x)[3])
end

function defaults(x::AbstractCompartmentSystem)
    build_toplevel(x)[4]
end

function get_systems(x::AbstractCompartmentSystem; rebuild = false)
    # collect channels + synapses + input systems
    #if rebuild || isempty(getfield(x, :systems))
        empty!(getfield(x, :systems))
        union!(getfield(x, :systems), getfield(x, :chans), getfield(x, :synapses), first.(getfield(x, :axial_conductance)))
    #end
    return getfield(x, :systems)
end

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

function Base.convert(::Type{ODESystem}, compartment::CompartmentSystem)
    required_states = Set{Num}()
    eqs, dvs, ps, defs, currents, syncurrents = build_toplevel(compartment)
    syss = convert.(ODESystem, collect(get_systems(compartment)))
    
    # Resolve in/out: "connect" / auto forward cell states to channels
    for x in get_channels(compartment) #union(get_channels(compartment), get_synapses(compartment))
        for inp in get_inputs(x)
            # gathering inputs
            # FIXME: add metadata indicating the level we expect to resolve at
            # e.g. Vₘ used by a synapse is resolved at the network level (an external source)
            # therefore it should remain unresolved when building the _host compartment_
            push!(required_states, inp)
            # all inputs are outright connected to the compartment top-level
            subinp = getproperty(x, tosymbol(inp, escape=false))
            push!(eqs, inp ~ subinp)
        end
    end
 
    int_modified = Set{Num}()
    foreach(x -> get_variables!(int_modified, x.lhs),  eqs)
   
    # TODO: filter required states from stimuli and extensions
    union!(required_states, setdiff(dvs, int_modified))
    # Resolve unavailable states
    for s in required_states
        # resolvedby(s) !== compartment && continue
        # Handled based on metadata of each state (for now just one case)
        if iscurrent(s) && isaggregate(s)
            push!(eqs, s ~ sum(filter(x -> getion(x) == getion(s), currents)))
            push!(dvs, s)
        end
        # ... other switch cases for required state handlers ...
    end
    return ODESystem(eqs, t, dvs, ps; systems = syss, defaults = defs, name = nameof(compartment))
end

#=
@enum StimulusTrigger ContinuousTrigger EdgeTriggered

abstract type Stimulus end

struct CurrentStimulus
    waveform
    offset::Current
    trigger::StimulusTrigger
    delay::TimeF64
    function CurrentStimulus(waveform; offset::Current = 0nA,
                             trigger::StimulusTrigger = ContinuousTrigger,
                             delay::Time = 0ms)
        return new(waveform, offset, trigger, delay)
    end
end

struct VoltageStimulus end
=#

