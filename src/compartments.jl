abstract type Geometry end
abstract type Sphere <: Geometry end
abstract type Cylinder <: Geometry end

# Extension building block
struct AuxConversion
    params::Vector{Num}
    eqs::Vector{Equation}
end

struct CompartmentSystem
    currents
    voltage
    chans::Vector{<:AbstractConductanceSystem}
    synapses::Vector{<:AbstractConductanceSystem}
    states::Vector
    params::Vector
    eqs
    defaults
    geometry
    systems
    observed
end

function CompartmentSystem(Vₘ, channels, aux, states, gradients, params; name, defaults) end

function CompartmentSystem(eqs, channels, states, params; gradients, name, aux)

    #Vₘ = MembranePotential()
    #params = @parameters cₘ aₘ
    
    grad_meta = getmetadata.(gradients, ConductorEquilibriumCtx)
    channel_systems = AbstractSystem[]
    eqs = Equation[] # equations must be a vector
    required_states = [] # states not produced or intrinsic (e.g. not currents or Vm)
    states = Any[Vₘ, Iapp, Isyn] # grow this as we discover/generate new states
    currents = []
    defaultmap = Pair[Iapp => ustrip(Float64, µA, holding),
                      aₘ => area,
                      Vₘ => ustrip(Float64, mV, V0),
                      Isyn => 0,
                      cₘ => ustrip(Float64, mF/cm^2, capacitance)]

    # By default, applied current is constant (just a bias/offset/holding current)
    # TODO: lift "stimulus" to a pass that happens at a higher level (i.e. after neuron
    # construction; with its own data type)
    if stimulus == nothing
        append!(eqs, [D(Iapp) ~ 0])
    else
        push!(eqs, Iapp ~ stimulus(t,Iapp))
    end

    # auxillary state transformations (e.g. net calcium current -> Ca concentration)
    if aux !== nothing
        for i in aux
            append!(params, i.params)
            # TODO: Probably can eliminate this block; MTK should find it automatically
            #for x in i.params
            #    hasdefault(x) && push!(defaultmap, x => getdefault(x))
            #end

            # gather all unique variables)
            inpvars = value.(vcat((get_variables(x.rhs) for x in i.eqs)...))
            unique!(inpvars)
            filter!(x -> !isparameter(x), inpvars) # exclude parameters
            append!(required_states, inpvars)

            # isolate states produced
            outvars = vcat((get_variables(x.lhs) for x in i.eqs)...)
            append!(states, outvars)
            append!(eqs, i.eqs)
            for j in outvars
                # FIXME: consider more consistent use of default variable ctx
                isconcentration(j) && push!(defaultmap, j => ustrip(Float64, µM, getconcentration(j).val))
            end
        end
    end

    # parse and build channel equations
    for chan in channels
        ion = chan.conducts
        sys = chan.sys
        push!(channel_systems, sys)

        # auto forward cell states to channels
        for inp in chan.inputs
            push!(required_states, inp)
            subinp = getproperty(sys, tosymbol(inp, escape=false))
            push!(eqs, inp ~ subinp)
            # Workaround for: https://github.com/SciML/ModelingToolkit.jl/issues/1013
            push!(defaultmap, subinp => inp)
        end

        # write the current equation state
        I = MembraneCurrent{ion}(name = nameof(sys), aggregate = false)
        push!(states, I)
        push!(currents, I)

        # for now, take the first reversal potential with matching ion type
        idx = findfirst(x -> x.ion == ion, grad_meta)
        Erev = gradients[idx]
        eq = [I ~ aₘ * sys.g * (Vₘ - Erev)]
        rhs = grad_meta[idx].val

        # check to see if reversal potential already defined
        if any(isequal(Erev, x) for x in states)
            append!(eqs, eq)
        else
            if typeof(rhs) <: Voltage
                push!(defaultmap, Erev => ustrip(Float64, mV, rhs))
                push!(params, Erev)
                append!(eqs, eq)
            else # symbolic/dynamic reversal potentials
                push!(eq, Erev ~ rhs)
                push!(states, Erev)
                rhs_vars = get_variables(rhs)
                filter!(x -> !isequal(x, value(Erev)), rhs_vars)
                rhs_ps = filter(x -> isparameter(x), rhs_vars)
                append!(eqs, eq)
                append!(required_states, rhs_vars)
            end
        end
    end

    required_states = unique(value.(required_states))
    states = Any[unique(value.(states))...]
    filter!(x -> !any(isequal(y, x) for y in states), required_states)

    if !isempty(required_states)
        newstateeqs = Equation[]
        for s in required_states
            # Handled based on metadata of each state (for now just one)
            if ismembranecurrent(s) && isaggregator(s)
                push!(newstateeqs, s ~ sum(filter(x -> iontype(x) == iontype(s), currents)))
                push!(states, s)
            end
            # ... other required state handlers ...
        end
        append!(eqs, newstateeqs)
    end

    # propagate default parameter values to channel systems
    vm_eq = D(Vₘ) ~ (Iapp - (+(currents..., Isyn)))/(aₘ*cₘ)
    push!(eqs, vm_eq)
    compartment = ODESystem(eqs, t, states, params; defaults = defaultmap, name = name)
    comp_sys = ModelingToolkit.compose(compartment, channel_systems)
    return Compartment{Sphere}(capacitance, channels, states, params, comp_sys)
end

const Soma = Compartment{Sphere}

function Compartment{Cylinder}() end
# should also be able to parse "collections" of compartments" that have an adjacency list/matrix
