
abstract type Geometry end

struct Sphere <: Geometry
    radius
end

function Sphere(; radius)
    Sphere(radius)
end

struct Point <: Geometry end

struct Unitless <: Geometry
    value
end

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
radius(::Union{Point,Unitless}) = 0.0

area(x::Sphere) = ustrip(Float64, cm^2, 4*π*radius(x)^2)
function area(x::Cylinder)
    ustrip(Float64, cm^2, 2*π*radius(x)*(height(x) + (x.open_ends ? 0µm : radius(x))))
end
area(::Point) = 1.0
area(x::Unitless) = x.value

"""
$(TYPEDEF)

A neuronal compartment.

$(TYPEDFIELDS)
"""
struct CompartmentSystem <: AbstractCompartmentSystem
    "Independent variabe. Defaults to time, ``t``."
    iv::Num
    "Voltage potential."
    voltage::Num
    "Membrane capacitance."
    capacitance::Num
    "Morphological geometry of the compartment."
    geometry::Geometry
    "Ionic conductances."
    chans::Set{AbstractConductanceSystem}
    "Equilibrium potentials belonging to ionic membrane conductances."
    channel_reversals::Set{Num}
    "Synaptic conductances."
    synapses::Set{AbstractConductanceSystem}
    "Equilibrium potentials belonging to synaptic conductances."
    synaptic_reversals::Set{Num}
    "Axial (intercompartmental) conductances."
    axial_conductance::Set{Tuple{AbstractConductanceSystem,Num}}
    "Experimental stimuli (for example, current injection)."
    stimuli::Vector{Equation}
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    defaults::Dict
    name::Symbol
    eqs::Vector{Equation}
    systems::Vector{AbstractTimeDependentSystem}
    observed::Vector{Equation}
    "Refers to the parent system when the compartment is a subcompartment in a `MultiCompartmentSystem`."
    parent::Ref{AbstractCompartmentSystem}
    function CompartmentSystem(iv, voltage, capacitance, geometry, chans, channel_reversals,
                               synapses, synaptic_reversals, axial_conductance, stimuli,
                               extensions, defaults, name, eqs, systems, observed, parent;
                               checks = false)
        if checks
        # placeholder
        end
        new(iv, voltage, capacitance, geometry, chans, channel_reversals, synapses,
            synaptic_reversals, axial_conductance, stimuli, extensions, defaults, name, eqs,
            systems, observed, parent)
    end
end

const Compartment = CompartmentSystem

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
    parent = Ref{AbstractCompartmentSystem}()
    return CompartmentSystem(t, Vₘ, cₘ, geometry, Set(channels), Set(reversals), Set(),
                             Set(), Set(), stimuli, extensions, Dict(), name, eqs, systems,
                             observed, parent)
end

# AbstractSystem interface extensions
get_geometry(x::AbstractCompartmentSystem) = getfield(x, :geometry)
area(x::AbstractCompartmentSystem) = only(@parameters aₘ = area(get_geometry(x)))
capacitance(x::AbstractCompartmentSystem) = getfield(x, :capacitance)
get_output(x::AbstractCompartmentSystem) = getfield(x, :voltage)

get_extensions(x::AbstractCompartmentSystem) = getfield(x, :extensions)
get_channels(x::AbstractCompartmentSystem) = getfield(x, :chans)
get_synapses(x::AbstractCompartmentSystem) = getfield(x, :synapses)
get_stimuli(x::AbstractCompartmentSystem) = getfield(x, :stimuli)

hasparent(x::CompartmentSystem) = isassigned(getfield(x, :parent))
parent(x::CompartmentSystem) = getfield(x, :parent)[]

function setparent!(child::CompartmentSystem, parent::AbstractCompartmentSystem)
    ref = getfield(child, :parent)
    ref[] = parent
    return nothing
end

function get_reversals(x::AbstractCompartmentSystem)
    return getfield(x, :channel_reversals), getfield(x, :synaptic_reversals) 
end

get_synaptic_reversals(x::AbstractCompartmentSystem) = get_reversals(x)[2]
get_channel_reversals(x::AbstractCompartmentSystem) = get_reversals(x)[1]

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
    
    filter!(x -> !isequal(x, t), rev_eq_vars)
    foreach(x -> isparameter(x) && push!(ps, x), rev_eq_vars)

    # Gather channel current equations
    for chan in get_channels(comp_sys)
        ion = permeability(chan)
        Erev = only(filter(x -> isequal(getion(x), ion), chan_revs))
        I = IonCurrent(ion, name = Symbol("I", nameof(chan)))
        g = renamespace(chan, get_output(chan))
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
        push!(eqs, I ~ g*aₘ*(childvm - Vₘ))
        push!(dvs, I)
        push!(dvs, childvm)
        push!(currents, -I)
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

    # Voltage equation
    push!(eqs, D(get_output(comp_sys)) ~ -sum(union(currents,syncurrents))/(cₘ*aₘ))

    return dvs, ps, eqs, defs, currents, syncurrents
end

# collect eqs including from extension + currents + reversals + Vₘ
function MTK.get_eqs(x::AbstractCompartmentSystem; rebuild = false)
    empty!(getfield(x, :eqs))
    union!(getfield(x, :eqs), build_toplevel(x)[3])
    return getfield(x, :eqs)
end

function MTK.get_states(x::AbstractCompartmentSystem)
    collect(build_toplevel(x)[1])
end

MTK.has_ps(x::CompartmentSystem) = true

# collect parameters from extension + currents + capacitance + area + reversals
function MTK.get_ps(x::AbstractCompartmentSystem)
    collect(build_toplevel(x)[2])
end

function MTK.defaults(x::AbstractCompartmentSystem)
    build_toplevel(x)[4]
end

# collect channels + synapses + input systems
function MTK.get_systems(x::AbstractCompartmentSystem; rebuild = false)
    empty!(getfield(x, :systems))
    union!(getfield(x, :systems), getfield(x, :chans), getfield(x, :synapses),
           first.(getfield(x, :axial_conductance)))
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
    dvs, ps, eqs, defs, currents, syncurrents = build_toplevel(compartment)
    syss = convert.(ODESystem, collect(get_systems(compartment)))
    
    # Resolve in/out: "connect" / auto forward cell states to channels
    for x in union(get_channels(compartment), get_synapses(compartment))
        for inp in get_inputs(x)
            isextrinsic(inp) && continue
            push!(required_states, inp)
            subinp = getproperty(x, tosymbol(inp, escape=false))
            push!(eqs, inp ~ subinp)
        end
    end
 
    internally_modified = Set{Num}()
    foreach(x -> get_variables!(internally_modified, x.lhs),  eqs)
   
    # TODO: also parse stimuli and extensions for required states
    union!(required_states, setdiff(dvs, internally_modified))

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
    return ODESystem(eqs, t, dvs, ps; systems = syss, defaults = defs,
                     name = nameof(compartment))
end

