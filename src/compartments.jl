
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
    area = 2*π*radius*(h + open_ends ? 0 : radius)
    return Cylinder(radius, length, area, open_ends)
end

length(x::Geometry) = hasfield(x, :length) ? getfield(x, :length) : nothing
radius(x::Geometry) = getfield(x, :radius)
radius(::Point) = 0.0

area(x::Sphere) = 4*π*radius(x)^2
area(x::Cylinder) = 2*π*radius(x)*(h + x.open_ends ? 0 : radius(x))
area(::Point) = 1.0

struct CompartmentSystem
    voltage # membrane voltage
    applied_current
    chans::Vector{<:AbstractConductanceSystem} # conductance systems giving rise to currents
    synapses::Vector{<:AbstractConductanceSystem} # synaptic conductance systems giving rise to synaptic currents
    states::Vector # all states
    params::Vector # all params
    eqs # equations _other than_ the voltage equation
    defaults::Dict
    geometry::Geometry # compartment geometry metadata (shape, dimensions, etc)
    aux_systems::Vector{ODESystem}
end

function CompartmentSystem(
    Vₘ,
    channels,
    reversals;
    capacitance = 1µF/cm^2,
    geometry::Geometry = Point(),
    aux::Vector{ODESystem} = ODESystem[],
    defaults = Dict(),
    name::Symbol = Base.gensym("Compartment")
) 
    # TODO: use built-in MTK units
    @parameters cₘ = ustrip(Float64, mF/cm^2, capacitance) aₘ = ustrip(Float64, cm^2, area(geometry))

    return CompartmentSystem(Vₘ, channels, states, params, systems)
end

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

