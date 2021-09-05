
struct CompartmentSystem
    voltage # membrane voltage
    currents # set of all currents
    chans       ::Vector{<:AbstractConductanceSystem} # conductance systems giving rise to currents
    synapses    ::Vector{<:AbstractConductanceSystem} # synaptic conductance systems giving rise to synaptic currents
    states      ::Vector # all states
    params      ::Vector # all params
    eqs # equations _other than_ the voltage equation
    defaults    ::Dict
    geometry # compartment geometry metadata (shape, dimensions, etc)
    systems # subsystems (e.g. ConductanceSystem)
    observed # placeholder; not implemented for now
end

function CompartmentSystem(
    Vₘ,
    channels,
    gradients,
    capacitance = 1µF/cm^2;
    geometry,
    aux::Vector{ODESystem} = ODESystem[],
    defaults = Dict(),
    name::Symbol = Base.gensym("Compartment")
) 

    @parameters cₘ aₘ
    grad_meta = getreversal.(gradients)
    
    eqs = Equation[]
    required_states = Set{Num}()
    dvs = Set{Num}([Vₘ])
    currents = Set{Num}()

    push!(defaults, aₘ => area)
    push!(defaults, cₘ => ustrip(Float64, mF/cm^2, capacitance))

    # auxillary state transformations (e.g. net calcium current -> Ca concentration)
    # ODESystems that extend the CompartmentSystem, which are flattened first.
    for auxsys in aux
        union!(eqs, equations(auxsys))
        union!(params, parameters(auxsys))
        union!(dvs, states(auxsys))
        defaults = merge(MTK.defaults(auxsys), defaults)
        aux_outputs = Set{Num}()

        for eq in equations(auxsys)
            modified_states!(aux_outputs, eq)
        end
        
        aux_inputs = setdiff(states(auxsys), aux_outputs)
        union!(required_states, aux_inputs)
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

    return CompartmentSystem(channels, states, params, systems)
end

