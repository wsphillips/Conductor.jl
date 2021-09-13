
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

struct CompartmentSystem
    ivs::Num # usually just t
    ## Intrinsic properties
    voltage::Num # symbol that represents membrane voltage
    capacitance::Num # specific membrane capacitance
    geometry::Geometry # compartment geometry metadata (shape, dimensions, etc)
    ## Dynamics
    currents::Set{Num} # container for static/dynamic currents
    chans::Set{AbstractConductanceSystem} # conductance systems
    channel_reversals::Set{Num}
    synapses::Set{AbstractConductanceSystem} # synaptic conductance systems
    synaptic_reversals::Set{Num}
    defaults::Dict
    aux_systems::Vector{ODESystem}
end

function CompartmentSystem(
    Vₘ,
    channels::Set{AbstractConductanceSystem},
    reversals;
    capacitance = 1µF/cm^2,
    geometry::Geometry = Point(),
    aux::Vector{ODESystem} = ODESystem[],
    name::Symbol = Base.gensym("Compartment")
) 
    @parameters cₘ = ustrip(Float64, mF/cm^2, capacitance)
    foreach(x -> isreversal(x) || throw("Invalid Equilibrium Potential"), reversals)
    return CompartmentSystem(t, Vₘ, cₘ, geometry, Set(), channels, Set(reversals), Set(),
                             Set(), Set(), Dict(), aux)
end

# AbstractSystem interface extensions
geometry(x::AbstractCompartmentSystem) = getfield(x, :geometry)
area(x::AbstractCompartmentSystem) = only(@parameters aₘ = area(geometry(x)))
capacitance(x::AbstractCompartmentSystem) = getfield(x, :capacitance)
output(x::AbstractCompartmentSystem) = getfield(x, :voltage)
get_auxsystems(x::AbstractCompartmentSystem) = getfield(x, :aux_systems)

function get_reversals(x::AbstractCompartmentSystem)
    return getfield(x, :channel_reversals), getfield(x, :synaptic_reversals) 
end

function build_toplevel(comp_sys::CompartmentSystem)
    dvs = Set{Num}()
    ps = Set{Num}()
    eqs = Equation[]
    
    aₘ = area(comp_sys)
    cₘ = capactiance(comp_sys)
    Vₘ = output(comp_sys)

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
    for chan in channels
        ion = permeability(chan)
        Erev = chan_revs[findfirst(x -> isequal(getion(x), ion), chan_revs)]
        I = IonCurrent(ion, name = nameof(chan))
        g = output(chan)
        # TODO: checks for whether we need area scaling
        push!(eqs, I ~ g*aₘ*(Vₘ - Erev))
        push!(dvs, I)
    end

    # Gather synaptic current equations
    for synapse in synapses
        ion = permeability(synapse)
        Esyn = syn_revs[findfirst(x -> isequal(getion(x), ion), syn_revs)]
        I = IonCurrent(ion, name = nameof(chan))
        g = output(synapse)
        push!(eqs, I ~ g*(Vₘ - Esyn))
        push!(dvs, I)
    end

    # Gather aux equations
    for aux in aux_systems
        union!(eqs, equations(aux))
        union!(ps, parameters(aux))
        union!(dvs, states(aux))
    end
    
    # stimuli pushed as currents?

    # voltage equation
    # push!(eqs, D(output(comp_sys)) ~ sum(currents)/cₘ*aₘ)

    return eqs, dvs, ps
end

# collect _top level_ eqs including from aux + currents + reversals + Vₘ
function get_eqs(comp_sys::AbstractCompartmentSystem)
    build_top_level(comp_sys)[1]
end

function get_states(x::AbstractCompartmentSystem)
    build_top_level(comp_sys)[2]
end

function get_ps(x::AbstractCompartmentSystem)
    # collect _top level_ parameters from aux + currents + capacitance + area + reversals

end

function get_systems(x::AbstractCompartmentSystem)
    # collect channels + synapses conductance systems
end

# System conversions

function Base.convert(ODESystem, compartment::CompartmentSystem)

    grad_meta = getreversal.(gradients)
    
    eqs = Equation[]
    dvs = Set{Num}([Vₘ])
    channel_currents = Set()
    required_states = Set{Num}()

    # auxillary state transformations (e.g. net calcium current -> Ca concentration)
    for auxsys in aux
        union!(eqs, equations(auxsys))
        union!(params, parameters(auxsys))
        union!(dvs, states(auxsys))
        defaults = merge(MTK.defaults(auxsys), defaults)
       
        # Push this to the end...
        #aux_outputs = Set{Num}()
        #for eq in equations(auxsys)
        #    modified_states!(aux_outputs, eq)
        #end
        #aux_inputs = setdiff(states(auxsys), aux_outputs)
        #union!(required_states, aux_inputs)
    end

    # parse and build channel equations
    for chan in channels
        # "connect" / auto forward cell states to channels
        for inp in getinputs(chan)
            push!(required_states, inp)
            subinp = getproperty(sys, tosymbol(inp, escape=false))
            push!(eqs, inp ~ subinp)
        end

        # symbol for current produced by channel
        ion = getion(chan)
        I = IonCurrent(ion, name = nameof(chan))
        push!(dvs, I)
        push!(currents, I)
        
        # Match to an available equilibirum potential
        # KLUDGE: take the first with matching ion type
        idx = findfirst(x -> getion(x) == ion, grad_meta)
        Erev = gradients[idx]
        
        # Current equation
        push!(eqs, I ~ aₘ * output(chan) * (Vₘ - Erev))

        if isparameter(Erev)
            push!(ps, Erev)
        else # dynamic eq potentials
            push!(dvs, Erev)
            rhs = getdefault(Erev)
            get_variables!(required_states, rhs)
            hasdefault(Erev) && push!(eqs, Erev ~ getdefault(Erev))
        end
    end

    # Resolve unavailable states
    setdiff!(required_states, dvs, ps, currents, reversals)

    for s in required_states
        # Handled based on metadata of each state (for now just one case)
        if iscurrent(s) && isaggregate(s)
            push!(eqs, s ~ sum(filter(x -> getion(x) == getion(s), currents)))
            push!(dvs, s)
        end
        # ... other switch cases for required state handlers ...
    end
    
    push!(eqs, D(Vₘ) ~ (Iapp - (+(currents..., Isyn)))/(aₘ*cₘ))

    return ODESystem(eqs, t, states, params; defaults = defaults, name = nameof(compartment))
end

