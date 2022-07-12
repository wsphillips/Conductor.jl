
get_output(x::AbstractConductanceSystem) = getfield(x, :output)
subscriptions(x::AbstractConductanceSystem) = getfield(x, :subscriptions)
isaggregate(x::AbstractConductanceSystem) = getfield(x, :aggregate)

function subscribe!(chan::AbstractConductanceSystem, comp)
    union!(subscriptions(chan), comp)
end

"""
$(TYPEDEF)

A model of conductance.

$(FIELDS)
"""
struct ConductanceSystem <: AbstractConductanceSystem
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
    "Conductance, ``g``, of the system."
    output::Num
    "Maximum conductance, ``\\overline{g}``."
    gbar::Num
    "Permeability of the conductance." 
    ion::IonSpecies
    aggregate::Bool # temp until better solution
    "Gating variables."
    gate_vars::Vector{AbstractGatingVariable}
    "Extrinsic sources of state (e.g. presynaptic compartments)."
    subscriptions::Vector{AbstractCompartmentSystem}
    """
    Additional systems to extend dynamics. Extensions are composed with the parent system
    during conversion to `ODESystem`.
    """
    extensions::Vector{ODESystem}
    inputs::Vector{Num}
    function ConductanceSystem(eqs, iv, states, ps, observed, name, systems, defaults,
            output, gbar, ion, aggregate, gate_vars, subscriptions, extensions, inputs;
                               checks = false)
        if checks
        #placeholder
        end
        new(eqs, iv, states, ps, observed, name, systems, defaults, output, gbar, ion,
            aggregate, gate_vars, subscriptions, extensions, inputs)
    end
end

const Conductance = ConductanceSystem

permeability(x::ConductanceSystem) = getfield(x, :ion)
get_gbar(x::ConductanceSystem) = getfield(x, :gbar)
get_gate_vars(x::ConductanceSystem) = getfield(x, :gate_vars)
"""
    ConductanceSystem(g, ion, gate_vars; <keyword arguments>)

Main constructor for `ConductanceSystem`.

# Arguments
- `gbar::Num`: Maximum conductance, ``\\overline{g}``.
- `aggregate::Bool`: determines whether the Conductance model should aggregate extrinsic
  sources of state instead of integrating them independently. Defaults to `false`.
- `defaults::Dict`: Default values for states and parameters.
- `name::Symbol`: Name of the system.
"""
function ConductanceSystem(g::Num,
                           ion::IonSpecies,
                           gate_vars::Vector{<:AbstractGatingVariable};
                           gbar::Num,
                           aggregate = false,
                           subscriptions = Vector{AbstractCompartmentSystem}(),
                           extensions::Vector{ODESystem} = ODESystem[],
                           defaults = Dict(),
                           name::Symbol = Base.gensym("Conductance"))

    gbar = setmetadata(gbar, ConductorMaxConductance, true)

    # Fields that will be generated
    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    inputs = Set{Num}()
    dvs = Set{Num}()
    ps = Set{Num}()
   
    # incomplete initialization
    cond_sys = ConductanceSystem(eqs, t, Num[], Num[], observed, name, systems, defaults, g,
                                 gbar, ion, aggregate, gate_vars, subscriptions, extensions,
                                 Num[])

    gate_var_outputs = Set{Num}()
    embed_defaults = Dict()
   
    isparameter(gbar) && push!(ps, gbar)

    for var in gate_vars
        vareqs = get_eqs(var, cond_sys)
        o = output(var)
        push!(isparameter(o) ? ps : gate_var_outputs, o)
        foreach(x -> get_variables!(inputs, x), vareqs)
        union!(eqs, vareqs)
    end

    for sym in inputs
        isparameter(sym) && push!(ps, sym)
        hasdefault(sym) && push!(embed_defaults, sym => getdefault(sym))
    end
    
    if isempty(gate_vars)
        push!(eqs, g ~ gbar)
    else
        push!(eqs, g ~ gbar * prod(output(x)^exponent(x) for x in gate_vars))
    end

    # Remove parameters + generated states
    setdiff!(inputs, ps, gate_var_outputs)
    union!(dvs, inputs, Set(g), gate_var_outputs)
    merge!(embed_defaults, defaults) 

    cond_sys = ConductanceSystem(eqs, t, collect(dvs), collect(ps), observed, name, systems,
                                 embed_defaults, g, gbar, ion, aggregate, gate_vars,
                                 subscriptions, extensions, collect(inputs))
    return cond_sys
end

function ConductanceSystem(x::ConductanceSystem;
                           g = get_output(x),
                           ion = permeability(x),
                           gate_vars = get_gate_vars(x),
                           gbar = get_gbar(x),
                           aggregate = isaggregate(x),
                           subscriptions = subscriptions(x),
                           extensions = get_extensions(x), 
                           defaults = get_defaults(x),
                           name = nameof(x))

    ConductanceSystem(g, ion, gate_vars;
                      gbar = gbar,
                      aggregate = aggregate,
                      subscriptions = subscriptions,
                      extensions = extensions,
                      defaults = defaults,
                      name = name)
end

get_inputs(x::ConductanceSystem) = getfield(x, :inputs)
get_extensions(x::AbstractConductanceSystem) = getfield(x, :extensions)

function Base.convert(::Type{ODESystem}, condsys::ConductanceSystem)
    dvs = states(condsys)
    ps  = parameters(condsys)
    eqs = equations(condsys)
    defs = get_defaults(condsys)
    sys = ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(condsys))
    #return extend(sys, get_extensions(sys))
    return sys
end

import ModelingToolkit: _eq_unordered

function Base.:(==)(sys1::ConductanceSystem, sys2::ConductanceSystem)
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

"""
    IonChannel(ion, gate_vars; <keyword arguments>)

An ionic membrane conductance.

# Arguments
- `max_g`: Default value for maximum conductance, ``\\overline{g}``.
- `extensions::Vector{ODESystem}`: Additional systems to extend dynamics. Extensions are
  composed with the parent system during conversion to `ODESystem`.
- `defaults::Dict`: Default values for states and parameters.
- `name::Symbol`: Name of the system.
"""
function IonChannel(ion::IonSpecies,
                    gate_vars::Vector{<:AbstractGatingVariable} = AbstractGatingVariable[];
                    max_g::Union{Num, SpecificConductance} = 0mS/cm^2,
                    extensions::Vector{ODESystem} = ODESystem[],
                    name::Symbol = Base.gensym("IonChannel"),
                    defaults = Dict())

    if max_g isa SpecificConductance
        gbar_val = ustrip(Float64, mS/cm^2, max_g)
        @parameters gbar
        push!(defaults, gbar => gbar_val)
    else
        gbar = max_g
        if hasdefault(gbar)
            gbar_val = getdefault(gbar)
            if gbar_val isa SpecificConductance
                gbar_val = ustrip(Float64, mS/cm^2, gbar_val)
                gbar = setdefault(gbar, gbar_val)
            end
        end
    end

    @variables g(t)
    g = setmetadata(g, ConductorUnits, mS/cm^2)
    gbar = setmetadata(gbar, ConductorUnits, mS/cm^2)

    ConductanceSystem(g, ion, gate_vars; gbar = gbar, name = name, defaults = defaults, 
                      extensions = extensions)
end

"""
    AxialConductance(gate_vars; <keyword arguments>)

A non-specific conductance between morphologically contiguous compartments.

# Arguments
- `max_g`: Maximum conductance, ``\\overline{g}``.
- `extensions::Vector{ODESystem}`: Additional systems to extend dynamics. Extensions are
  composed with the parent system during conversion to `ODESystem`.
- `defaults::Dict`: Default values for states and parameters.
- `name::Symbol`: Name of the system.
"""
function AxialConductance(gate_vars::Vector{<:AbstractGatingVariable} = AbstractGatingVariable[];
                          max_g = 0mS/cm^2, extensions::Vector{ODESystem} = ODESystem[],
                          name::Symbol = Base.gensym("Axial"), defaults = Dict())
    
    IonChannel(Leak, gate_vars, max_g = max_g, extensions = extensions, name = name,
               defaults = defaults)
end

"""
    SynapticChannel(ion, gate_vars; <keyword arguments>)

A synaptically activated conductance. Depends on extrinsic (i.e. presynaptic) state.

# Arguments
- `max_g`: Maximum conductance, ``\\overline{g}``.
- `extensions::Vector{ODESystem}`: Additional systems to extend dynamics. Extensions are
  composed with the parent system during conversion to `ODESystem`.
- `aggregate::Bool`: whether the Conductance model should aggregate extrinsic sources of
  state instead of integrating them independently. Defaults to `false`.
- `defaults::Dict`: Default values for states and parameters.
- `name::Symbol`: Name of the system.
"""
function SynapticChannel(ion::IonSpecies,
                         gate_vars::Vector{<:AbstractGatingVariable} = AbstractGatingVariable[];
                         max_s::Union{Num, ElectricalConductance} = 0mS,
                         extensions::Vector{ODESystem} = ODESystem[], aggregate = false,
                         name::Symbol = Base.gensym("SynapticChannel"), defaults = Dict())

    if max_s isa ElectricalConductance
        sbar_val = ustrip(Float64, mS, max_s)
        @parameters sbar = sbar_val
        push!(defaults, sbar => sbar_val)
    else
        sbar = max_s
        if hasdefault(sbar)
            sbar_val = getdefault(sbar)
            if sbar_val isa ElectricalConductance
                sbar_val = ustrip(Float64, mS, sbar_val)
                sbar = setdefault(sbar, sbar_val)
            end
        end
    end

    @variables s(t)
    s = setmetadata(s, ConductorUnits, mS)
    ConductanceSystem(s, ion, gate_vars; gbar = sbar, name = name, defaults = defaults,
                      aggregate = aggregate, extensions = extensions)
end

function (cond::AbstractConductanceSystem)(newgbar::Num)
    hasdefault(newgbar) || throw("Symbolic has no default value")
    gbar_val = getdefault(newgbar)
    return cond(gbar_val)
end

function (cond::AbstractConductanceSystem)(newgbar::Quantity)
    g = get_output(cond)
    outunits = getmetadata(g, ConductorUnits)
    if dimension(outunits) !== dimension(newgbar)
        @error "Input Dimensions do not match output of ConductanceSystem"
    end
    gbar_val = ustrip(Float64, outunits, newgbar)
    return cond(gbar_val)
end

function (cond::AbstractConductanceSystem)(gbar_val::Real)
    gbar_sym = setdefault(get_gbar(cond), gbar_val)
    newcond = @set cond.defaults = Dict(get_defaults(cond)..., gbar_sym => gbar_val)
    return newcond
end

# Old macros -- needs updating
#= 
macro ionchannel(chan::Expr, ex::Expr, gbar)
    gbar.args[1] != :gbar && throw("Please use `gbar` to define a channel conductance.")
    name, ion, gates = _make_ionchannel(chan, ex)
    IonChannel(eval(ion), gates; max_g = eval(gbar.args[2]), name = name)
end

macro ionchannel(chan::Expr, ex::Expr)
    name, ion, gates = _make_ionchannel(chan, ex)
    IonChannel(eval(ion), gates; name = name)
end

function _make_ionchannel(chan::Expr, ex::Expr)
    ex = MacroTools.striplines(ex)
    !@capture(chan, name_Symbol{ion_}) && throw("An ion type must be given `name{I<:Ion}`")
    gates = Gate[]
    for gate in ex.args
      push!(gates, eval(gate))
    end
    name, ion, gates
end
=#
