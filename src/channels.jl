# Conductance types (conductance as in "g")
abstract type AbstractConductanceSystem <: AbstractTimeDependentSystem end

# Linear Ohmic/Nernst vs non-linear GHK
@enum IVCurvature Linear Rectifying

permeability(x::AbstractConductanceSystem) = getfield(x, :ion)
get_inputs(x::AbstractConductanceSystem) = getfield(x, :inputs)
get_output(x::AbstractConductanceSystem) = getfield(x, :output)
# Abstract types without parametrics
struct ConductanceSystem{S<:AbstractTimeDependentSystem} <: AbstractConductanceSystem
    
    iv
    output::Num # 'g' by default 
    ion::IonSpecies # ion permeability
    gate_vars::Vector{AbstractGatingVariable}
    inputs::Set{Num} # required inputs
    sys::S
    linearity::IVCurvature
    transmatrix::Union{Matrix, Nothing}
    name::Symbol
    eqs::Vector{Equation}
    systems::Vector{AbstractTimeDependentSystem}
    function ConductanceSystem(iv, output, ion, gate_vars, inputs, sys, linearity, transmatrix,
                               name, eqs, systems; checks = false)
        if checks
        #placeholder
        end
        new{typeof(sys)}(iv, output, ion, gate_vars, inputs, sys, linearity, transmatrix, name, eqs, systems)
    end
end

const Conductance = ConductanceSystem

import ModelingToolkit: _eq_unordered

# This should trigger a rebuild instead
Base.convert(::Type{ODESystem}, x::ConductanceSystem{ODESystem}) = getfield(x, :sys)

function ModelingToolkit.rename(x::ConductanceSystem, name)
    xcopy = deepcopy(x)
    @set! xcopy.sys.name = name
    @set xcopy.name = name
end

# Forward getters to internal system
for prop in [
             :eqs
             :noiseeqs
             :iv
             :states
             :ps
             :var_to_name
             :ctrls
             :defaults
             :observed
             :tgrad
             :jac
             :ctrl_jac
             :Wfact
             :Wfact_t
             :systems
             :structure
             :op
             :equality_constraints
             :inequality_constraints
             :controls
             :loss
             :bcs
             :domain
             :ivs
             :dvs
             :connector_type
             :connections
             :preface
             :torn_matching
             :tearing_state
             :substitutions
            ]
    fname1 = Symbol(:get_, prop)
    fname2 = Symbol(:has_, prop)
    @eval begin
        $fname1(x::ConductanceSystem) = getfield(getfield(x, :sys), $(QuoteNode(prop)))
        $fname2(x::ConductanceSystem) = isdefined(getfield(x, :sys), $(QuoteNode(prop)))
    end
end

MTK.getvar(x::ConductanceSystem, name::Symbol; namespace::Bool) = MTK.getvar(getfield(x, :sys), name, namespace = namespace)

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

function build_gate_eq(var::Gate{<:Union{AlphaBeta,SteadyStateTau}})
    x, x∞, τₓ = output(var), steadystate(var), timeconstant(var)
    return D(x) ~ inv(τₓ)*(x∞ - x)
end

build_gate_eq(var::Gate{<:Union{SteadyState,ConstantValue}}) = output(var) ~ steadystate(var)

function ConductanceSystem(g::Num, ion::IonSpecies, gate_vars::Vector{<:AbstractGatingVariable};
        gbar::Num, linearity::IVCurvature = Linear, extensions::Vector{ODESystem} = ODESystem[],
                           defaults = Dict(), name::Symbol = Base.gensym("Conductance"))

    eqs = Equation[]
    inputs = Set{Num}()
    gate_var_outputs = Set{Num}()
    embed_defaults = Dict()
    params = Set{Num}()
    
    gbar = setmetadata(gbar, ConductorMaxConductance, true)
    isparameter(gbar) && push!(params, gbar)

    for var in gate_vars
        eq = build_gate_eq(var)
        o = output(var)
        push!(isparameter(o) ? params : gate_var_outputs, o)
        get_variables!(inputs, eq)
        push!(eqs, eq)
    end

    for sym in inputs
        isparameter(sym) && push!(params, sym)
        hasdefault(sym) && push!(embed_defaults, sym => getdefault(sym))
    end
    
    # Remove parameters + generated states
    setdiff!(inputs, params, gate_var_outputs)
    if length(gate_vars) > 0
        push!(eqs, g ~ gbar * prod(hasexponent(x) ? output(x)^exponent(x) : output(x) for x in gate_vars))
    else
        push!(eqs, g ~ gbar)
    end
    sys = ODESystem(eqs, t, union(inputs, g, gate_var_outputs), params;
                    defaults = merge(embed_defaults, defaults), name = name)
    
    for ext in extensions
        sys = extend(sys, ext)
    end

    return ConductanceSystem(t, g, ion, gate_vars, inputs, sys, linearity, nothing, name, eqs, Vector{ODESystem}[])
end

function IonChannel(ion::IonSpecies,
        gate_vars::Vector{<:AbstractGatingVariable} = AbstractGatingVariable[];
                    max_g::Union{Num, SpecificConductance} = 0mS/cm^2,
                    extensions::Vector{ODESystem} = ODESystem[],
                    name::Symbol = Base.gensym("IonChannel"),
                    linearity::IVCurvature = Linear, defaults = Dict())
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
                      extensions = extensions, linearity = linearity)
end

function SynapticChannel(ion::IonSpecies,
        gate_vars::Vector{<:AbstractGatingVariable} = AbstractGatingVariable[];
                         max_s::Union{Num, ElectricalConductance} = 0mS,
                         extensions::Vector{ODESystem} = ODESystem[],
                         name::Symbol = Base.gensym("SynapticChannel"),
                         linearity::IVCurvature = Linear, defaults = Dict())
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
                      gbar = sbar, name = name, defaults = defaults,
                      extensions = extensions, linearity = linearity)
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
        if getmetadata(param, ConductorMaxConductance, false)
            gbar_sym = param
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
