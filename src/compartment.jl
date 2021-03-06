
"""
$(TYPEDEF)

A neuronal compartment.

$(TYPEDFIELDS)
"""
struct CompartmentSystem <: AbstractCompartmentSystem
    # MTK fields
    eqs::Vector{Equation}
    "Independent variabe. Defaults to time, ``t``."
    iv::Num
    states::Vector{Num}
    ps::Vector{Num}
    observed::Vector{Equation}
    name::Symbol
    systems::Vector{AbstractTimeDependentSystem}
    defaults::Dict
    # Conductor fields
    "Voltage potential."
    voltage::Num
    "Membrane capacitance."
    capacitance::Num
    "Morphological geometry of the compartment."
    geometry::Geometry
    "Ionic conductances."
    chans::Vector{AbstractConductanceSystem}
    "Equilibrium potentials belonging to ionic membrane conductances."
    channel_reversals::Vector{Num}
    "Synaptic conductances."
    synapses::Vector{AbstractConductanceSystem}
    "Equilibrium potentials belonging to synaptic conductances."
    synaptic_reversals::Vector{Num}
    "Axial (intercompartmental) conductances."
    axial_conductance::Vector{Tuple{AbstractConductanceSystem,Num}}
    "Experimental stimuli (for example, current injection)."
    stimuli::Vector{Equation}
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    "Refers to the parent system when the compartment is a subcompartment in a `MultiCompartmentSystem`."
    parent::Ref{AbstractCompartmentSystem}
    function CompartmentSystem(eqs, iv, states, ps, observed, name, systems, defaults,
                               voltage, capacitance, geometry, chans, channel_reversals,
                               synapses, synaptic_reversals, axial_conductance, stimuli,
                               extensions, parent;
                               checks = false)
        if checks
        # placeholder
        foreach(x -> isreversal(x) || throw("Invalid Equilibrium Potential"), channel_reversals)
        foreach(x -> iscurrent(x.lhs) || throw("Invalid current stimulus"), stimuli)

        end
        new(eqs, iv, states, ps, observed, name, systems, defaults,
                               voltage, capacitance, geometry, chans, channel_reversals,
                               synapses, synaptic_reversals, axial_conductance, stimuli,
                               extensions, parent)
    end
end

const Compartment = CompartmentSystem

"""
    CompartmentSystem(V???, channels, reversals; <keyword arguments>)

# Arguments
- `capacitance::SpecificCapacitance`: The capacitance of the compartment given in Farads 
  per unit area (e.g. ??F/cm^2).
- `geometry::Geometry`: Morphological geometry of the compartment.
- `extensions::Vector{ODESystem}`: Additional systems to extend dynamics. Extensions are
  composed with the parent system during conversion to `ODESystem`.
- `stimuli::Vector{Equation}`: 
- `name::Symbol`: Name of the system.
"""
function CompartmentSystem(
    V???::Num,
    channels,
    channel_reversals;
    capacitance = 1??F/cm^2,
    geometry::Geometry = Point(),
    extensions::Vector{ODESystem} = ODESystem[],
    stimuli::Vector{Equation} = Equation[],
    defaults = Dict(),
    name::Symbol = Base.gensym("Compartment")
) 
    
    @parameters c??? = ustrip(Float64, mF/cm^2, capacitance)
    Set(channels)
    Set(channel_reversals)
    synaptic_channels = Set{AbstractConductanceSystem}()
    synaptic_reversals = Set{Num}()
    axial_conductance = Set{Tuple{AbstractConductanceSystem,Num}}()
    parent = Ref{AbstractCompartmentSystem}()

    return CompartmentSystem(V???, c???, geometry, channels, channel_reversals, 
                             synaptic_channels, synaptic_reversals, axial_conductance, 
                             stimuli, extensions, parent, name, defaults)
end

# takes all user specified fields and generates remaining fields
function CompartmentSystem(V???, c???, geometry,
                           channels, channel_reversals,
                           synaptic_channels, synaptic_reversals,
                           axial_conductance,
                           stimuli, extensions, parent, name, defaults)

    # generated
    eqs = Equation[]
    systems = union(channels, synaptic_channels, first.(axial_conductance))
    observed = Equation[]
    dvs = Set{Num}()
    ps  = Set{Num}()
    currents = Set()
    defs = Dict()

    @parameters a??? = area(geometry)
    push!(dvs, V???)
    push!(ps, a???)
    push!(ps, c???)

    # Check for possible dynamic reversals
    reversal_equation_vars = Set{Num}()
    for Erev in union(channel_reversals, synaptic_reversals, last.(axial_conductance))
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
    
    paired_channels = zip(channels,
                          broadcast(chan -> find_reversal(chan, channel_reversals), channels))
    paired_synapses = zip(synaptic_channels,
                          broadcast(chan -> find_reversal(chan, synaptic_reversals), synaptic_channels))
    paired_conductances = union(axial_conductance, paired_channels, paired_synapses)

    for (chan, Erev) in paired_conductances
        I = IonCurrent(chan)
        g = renamespace(chan, get_output(chan))
        push!(eqs, I ~ g*(V??? - Erev)*(1*getmetadata(g, ConductorUnits) isa SpecificConductance ? a??? : 1))
        push!(dvs, I)
        push!(currents, I)
    end

    for stimulus in stimuli
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
    push!(eqs, D(V???) ~ -sum(currents)/(c???*a???))

    # compose extensions
    for extension in extensions
        union!(eqs, equations(extension))
        union!(ps, parameters(extension))
        union!(dvs, states(extension))
    end

### previously part of convert(ODESystem, x)
#
    required_states = Set{Num}()
    # Resolve in/out: "connect" / auto forward cell states to channels
    for x in union(channels, synaptic_channels)
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
    component_currents = filter(x -> iscurrent(x) && !isaggregate(x), union(dvs,ps))
    # Resolve unavailable states
    for s in required_states
        # resolvedby(s) !== compartment && continue
        # Handled based on metadata of each state (for now just one case)
        if iscurrent(s) && isaggregate(s)
            push!(eqs, s ~ sum(filter(x -> getion(x) == getion(s), component_currents)))
            push!(dvs, s)
        end
        # ... other switch cases for required state handlers ...
    end
   
    merge!(defs, defaults)

    return  CompartmentSystem(eqs, t, collect(dvs), collect(ps), observed, name,
                              collect(systems), defs, V???, c???, geometry, collect(channels),
                              collect(channel_reversals), collect(synaptic_channels),
                              collect(synaptic_reversals), collect(axial_conductance),
                              stimuli, extensions, parent)
end


function CompartmentSystem(sys::CompartmentSystem;
                           V??? = get_output(sys),
                           c??? = capacitance(sys),
                           geometry = get_geometry(sys),
                           channels = get_channels(sys),
                           channel_reversals = get_channel_reversals(sys),
                           synaptic_channels = get_synapses(sys), 
                           synaptic_reversals = get_synaptic_reversals(sys),
                           axial_conductance = get_axial_conductance(sys),
                           stimuli = get_stimuli(sys),
                           extensions = get_extensions(sys),
                           parent = parent(sys),
                           name = nameof(sys),
                           defaults = get_defaults(sys))
    CompartmentSystem(V???, c???, geometry,
                           channels, channel_reversals,
                           synaptic_channels, synaptic_reversals,
                           axial_conductance,
                           stimuli, extensions, parent, name, defaults)
end


# AbstractSystem interface extensions
get_geometry(x::AbstractCompartmentSystem) = getfield(x, :geometry)
area(x::AbstractCompartmentSystem) = only(@parameters a??? = area(get_geometry(x)))
capacitance(x::AbstractCompartmentSystem) = getfield(x, :capacitance)
get_output(x::AbstractCompartmentSystem) = getfield(x, :voltage)

get_extensions(x::AbstractCompartmentSystem) = getfield(x, :extensions)
get_channels(x::AbstractCompartmentSystem) = getfield(x, :chans)
get_synapses(x::AbstractCompartmentSystem) = getfield(x, :synapses)
get_stimuli(x::AbstractCompartmentSystem) = getfield(x, :stimuli)

hasparent(x::CompartmentSystem) = isassigned(getfield(x, :parent))
parent(x::CompartmentSystem) = getfield(x, :parent)

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

