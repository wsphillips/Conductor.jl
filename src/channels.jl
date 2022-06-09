
# Linear Ohmic/Nernst vs non-linear GHK
# @enum IVCurvature Linear Rectifying

get_output(x::AbstractConductanceSystem) = getfield(x, :output)
subscriptions(x::AbstractConductanceSystem) = getfield(x, :subscriptions)
isaggregate(x::AbstractConductanceSystem) = getfield(x, :aggregate)

function subscribe!(chan::AbstractConductanceSystem, comp::AbstractCompartmentSystem)
    push!(subscriptions(chan), comp)
end

struct ConductanceSystem <: AbstractConductanceSystem
    iv::Num
    output::Num # 'g' by default 
    gbar::Num # max conductance parameter
    ion::IonSpecies # ion permeability
    aggregate::Bool # temp until better solution
    gate_vars::Vector{AbstractGatingVariable}
    subscriptions::Set{AbstractCompartmentSystem}
    extensions::Vector{ODESystem}
    defaults::Dict
    name::Symbol
    eqs::Vector{Equation}
    systems::Vector{AbstractTimeDependentSystem}
    observed::Vector{Equation}
    function ConductanceSystem(iv, output, gbar, ion, aggregate, gate_vars, subscriptions, extensions,
                               defaults, name, eqs, systems, observed; checks = false)
        if checks
        #placeholder
        end
        new(iv, output, gbar, ion, aggregate, gate_vars, subscriptions, extensions, defaults, name,
            eqs, systems, observed)
    end
end

const Conductance = ConductanceSystem

permeability(x::ConductanceSystem) = getfield(x, :ion)

function ConductanceSystem(g::Num, ion::IonSpecies,
                           gate_vars::Vector{<:AbstractGatingVariable}; gbar::Num, aggregate = false,
                           extensions::Vector{ODESystem} = ODESystem[], defaults = Dict(),
                           name::Symbol = Base.gensym("Conductance"))

    gbar = setmetadata(gbar, ConductorMaxConductance, true)
    eqs = Equation[]
    systems = AbstractTimeDependentSystem[]
    observed = Equation[]
    subscriptions = Set{AbstractCompartmentSystem}()

    return ConductanceSystem(t, g, gbar, ion, aggregate, gate_vars, subscriptions, extensions,
                             defaults, name, eqs, systems, observed)
end

function build_toplevel!(dvs, ps, eqs, defs, comp_sys::ConductanceSystem)
    
    inputs = Set{Num}()
    gate_var_outputs = Set{Num}()
    embed_defaults = Dict()
    
    gbar = getfield(comp_sys, :gbar) # needs getter?
    isparameter(gbar) && push!(ps, gbar)
    g = get_output(comp_sys) 
    gate_vars = getfield(comp_sys, :gate_vars) # needs getter?

    for var in gate_vars
        vareqs = get_eqs(var, comp_sys)
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

    merge!(defs, merge(embed_defaults, getfield(comp_sys, :defaults))) 

    return dvs, ps, eqs, defs, inputs
end

get_inputs(x::ConductanceSystem) = build_toplevel(x)[5]

function get_eqs(x::AbstractConductanceSystem)
    empty!(getfield(x, :eqs))
    union!(getfield(x, :eqs), build_toplevel(x)[3])
    return getfield(x, :eqs)
end

function get_states(x::AbstractConductanceSystem)
    collect(build_toplevel(x)[1])
end

MTK.has_ps(x::ConductanceSystem) = true

function get_ps(x::AbstractConductanceSystem)
    collect(build_toplevel(x)[2])
end

function defaults(x::AbstractConductanceSystem)
    build_toplevel(x)[4]
end

function get_systems(x::AbstractConductanceSystem)
    #empty!(getfield(x, :systems))
    #union!(getfield(x, :systems), getfield(x, :chans), getfield(x, :synapses), first.(getfield(x, :axial_conductance)))
    return getfield(x, :systems)
end

get_extensions(x::AbstractConductanceSystem) = getfield(x, :extensions)

function Base.convert(::Type{ODESystem}, condsys::ConductanceSystem)
    dvs, ps, eqs, defs, _ = build_toplevel(condsys)
    sys = ODESystem(eqs, t, dvs, ps; defaults = defs, name = nameof(condsys))
    #return extend(sys, get_extensions(sys))
    return sys
end

#function ModelingToolkit.rename(x::ConductanceSystem, name)
#    xcopy = deepcopy(x)
#    @set! xcopy.sys.name = name
#    @set xcopy.name = name
#end

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
    g = setmetadata(g, ConductorUnits, mS/cm^2) # TODO: rework with MTK native unit system
    ConductanceSystem(g, ion, gate_vars;
                      gbar = gbar, name = name, defaults = defaults, 
                      extensions = extensions)
end

function AxialConductance(gate_vars::Vector{<:AbstractGatingVariable} = AbstractGatingVariable[];
                          max_g = 0mS/cm^2, extensions::Vector{ODESystem} = ODESystem[],
                          name::Symbol = Base.gensym("Axial"), defaults = Dict())
    
    IonChannel(Leak, gate_vars, max_g = max_g, extensions = extensions, name = name,
               defaults = defaults)
end

function SynapticChannel(ion::IonSpecies,
                         gate_vars::Vector{<:AbstractGatingVariable} = AbstractGatingVariable[];
                         max_s::Union{Num, ElectricalConductance} = 0mS,
                         extensions::Vector{ODESystem} = ODESystem[], aggregate = false,
                         name::Symbol = Base.gensym("SynapticChannel"),
                         defaults = Dict())
    # to make generic, check for <:Quantity then write a
    # unit specific "strip" method 
    if max_s isa ElectricalConductance
        sbar_val = ustrip(Float64, mS, max_s)
        @parameters sbar
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
    s = setmetadata(s, ConductorUnits, mS) # TODO: rework iwth MTK native unit system
    ConductanceSystem(s, ion, gate_vars;
                      gbar = sbar, name = name, defaults = defaults, aggregate = aggregate,
                      extensions = extensions)
end

function (cond::AbstractConductanceSystem)(newgbar::Quantity)
    
    newcond = deepcopy(cond)
    g = get_output(newcond)
    outunits = getmetadata(g, ConductorUnits)

    if dimension(outunits) !== dimension(newgbar)
        @error "Input Dimensions do not match output of ConductanceSystem"
    end

    gbar_val = ustrip(Float64, outunits, newgbar)
    local gbar_sym::Num

    for param in parameters(cond)
        if getmetadata(unwrap(param), ConductorMaxConductance, false)
            gbar_sym = wrap(param)
            break
        end
    end

    push!(get_defaults(newcond), gbar_sym => gbar_val)

    return newcond
end

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
