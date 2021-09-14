
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
    length
    open_ends::Bool
end

function Cylinder(; radius, length, open_ends = true)
    return Cylinder(radius, length, open_ends)
end

length(x::Geometry) = hasfield(x, :length) ? getfield(x, :length) : nothing
radius(x::Geometry) = getfield(x, :radius)
radius(::Point) = 0.0

area(x::Sphere) = 4*π*radius(x)^2
area(x::Cylinder) = 2*π*radius(x)*(h + x.open_ends ? 0 : radius(x))
area(::Point) = 1.0

# stub -- TODO: implement datastructure for describing electrode stimuli
struct Stimulus end

struct CompartmentSystem
    ivs::Num # usually just t
    ## Intrinsic properties
    voltage::Num # symbol that represents membrane voltage
    capacitance::Num # specific membrane capacitance
    geometry::Geometry # compartment geometry metadata (shape, dimensions, etc)
    ## Dynamics
    chans::Set{AbstractConductanceSystem} # conductance systems
    channel_reversals::Set{Num}
    synapses::Set{AbstractConductanceSystem} # synaptic conductance systems
    synaptic_reversals::Set{Num}
    stimuli::Vector{Stimulus}
    extensions::Vector{ODESystem}
    defaults::Dict
end

function CompartmentSystem(
    Vₘ,
    channels::Set{AbstractConductanceSystem},
    reversals;
    capacitance = 1µF/cm^2,
    geometry::Geometry = Point(),
    extensions::Vector{ODESystem} = ODESystem[],
    name::Symbol = Base.gensym("Compartment")
) 
    @parameters cₘ = ustrip(Float64, mF/cm^2, capacitance)
    foreach(x -> isreversal(x) || throw("Invalid Equilibrium Potential"), reversals)
    return CompartmentSystem(t, Vₘ, cₘ, geometry, Set(), channels, Set(reversals), Set(),
                             Set(), Set(), Dict(), extensions)
end

# AbstractSystem interface extensions
get_geometry(x::AbstractCompartmentSystem) = getfield(x, :geometry)
area(x::AbstractCompartmentSystem) = only(@parameters aₘ = area(get_geometry(x)))
capacitance(x::AbstractCompartmentSystem) = getfield(x, :capacitance)
output(x::AbstractCompartmentSystem) = getfield(x, :voltage)

get_extensions(x::AbstractCompartmentSystem) = getfield(x, :extensions)
get_channels(x::AbstractCompartmentSystem) = getfield(x, :chans)
get_synapses(x::AbstractCompartmentSystem) = getfield(x, :synapses)

function get_inputs(x::AbstractCompartmentSystem) end 
function inputs(x::AbstractCompartmentSystem) end

function get_reversals(x::AbstractCompartmentSystem)
    return getfield(x, :channel_reversals), getfield(x, :synaptic_reversals) 
end

function build_toplevel(comp_sys)
    dvs = Set{Num}()
    ps = Set{Num}()
    eqs = Equation[]
    defs = Dict()
    build_toplevel!(dvs, ps, eqs, defs, comp_sys)
    return eqs, dvs, ps, defs
end

function build_toplevel!(dvs, ps, eqs, defs, comp_sys::CompartmentSystem)
    aₘ = area(comp_sys)
    cₘ = capactiance(comp_sys)
    Vₘ = output(comp_sys)

    currents = Set{Num}()

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
    for chan in get_channels(x)
        ion = permeability(chan)
        Erev = chan_revs[findfirst(x -> isequal(getion(x), ion), chan_revs)]
        I = IonCurrent(ion, name = nameof(chan))
        g = output(chan)
        # TODO: checks for whether we need area scaling
        push!(eqs, I ~ g*aₘ*(Vₘ - Erev))
        push!(dvs, I)
        push!(currents, I)
        merge!(defs, defaults(chan))
    end

    # Gather synaptic current equations
    for synapse in get_synapses(x)
        ion = permeability(synapse)
        Esyn = syn_revs[findfirst(x -> isequal(getion(x), ion), syn_revs)]
        I = IonCurrent(ion, name = nameof(chan))
        g = output(synapse)
        push!(eqs, I ~ g*(Vₘ - Esyn))
        push!(dvs, I)
        push!(currents, I)
        merge!(defs, defaults(chan))
    end

    # Gather extension equations
    for extension in get_extensions(x)
        union!(eqs, equations(extension))
        union!(ps, parameters(extension))
        union!(dvs, states(extension))
        merge!(defs, defaults(extension))
    end
    
    # FIXME: separate field for applied current(s)?

    # voltage equation
    push!(eqs, D(output(comp_sys)) ~ -sum(currents)/cₘ*aₘ)
end

# collect _top level_ eqs including from extension + currents + reversals + Vₘ
function get_eqs(x::AbstractCompartmentSystem)
    build_top_level(x)[1]
end

function get_states(x::AbstractCompartmentSystem)
    build_top_level(x)[2]
end

function get_ps(x::AbstractCompartmentSystem)
    # collect _top level_ parameters from extension + currents + capacitance + area + reversals
    build_top_level(x)[3]
end

function defaults(x::AbstractCompartmentSystem)
    build_top_level(x)[4]
end

function get_systems(x::AbstractCompartmentSystem)
    # collect channels + synapses + input systems
    union(x.chans, x.synapses)
end

function Base.convert(ODESystem, compartment::CompartmentSystem)

    required_states = Set{Num}()
    eqs, dvs, ps, defs = build_toplevel(compartment)

    # Resolve in/out: "connect" / auto forward cell states to channels
    for x in union(get_channels(compartment), get_synapses(compartment))
        for inp in get_inputs(x)
            # gathering inputs
            push!(required_states, inp)
            # all inputs are outright connected to the compartment top-level
            subinp = getproperty(sys, tosymbol(inp, escape=false))
            push!(eqs, inp ~ subinp)
        end
    end
    
    setdiff!(required_states, dvs)

    # Resolve unavailable states
    for s in required_states
        # Handled based on metadata of each state (for now just one case)
        if iscurrent(s) && isaggregate(s)
            push!(eqs, s ~ sum(filter(x -> get_ion(x) == get_ion(s), currents)))
            push!(dvs, s)
        end
        # ... other switch cases for required state handlers ...
    end

    return ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(compartment))
end

