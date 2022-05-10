# Conductance types (conductance as in "g")
abstract type AbstractConductanceSystem <: AbstractTimeDependentSystem end

# Linear Ohmic/Nernst vs non-linear GHK
@enum IVCurvature Linear Rectifying

permeability(x::AbstractConductanceSystem) = getfield(x, :ion)
get_inputs(x::AbstractConductanceSystem) = getfield(x, :inputs)
get_output(x::AbstractConductanceSystem) = getfield(x, :output)
# Abstract types without parametrics
struct ConductanceSystem{S<:AbstractTimeDependentSystem} <: AbstractConductanceSystem
    output::Num # 'g' by default 
    ion::IonSpecies # ion permeability
    gate_vars::Vector{<:AbstractGatingVariable}
    inputs::Set{Num} # required inputs
    sys::S
    linearity::IVCurvature
    transmatrix::Union{Matrix, Nothing}
    name::Symbol
end

Base.convert(::Type{ODESystem}, x::ConductanceSystem{ODESystem}) = getfield(x, :sys)

# Forward getters to internal system
MTK.get_systems(x::AbstractConductanceSystem) = get_systems(getfield(x, :sys))
MTK.get_eqs(x::AbstractConductanceSystem) = get_eqs(getfield(x, :sys))
MTK.get_dvs(x::AbstractConductanceSystem) = get_dvs(getfield(x, :sys))
MTK.has_ps(x::AbstractConductanceSystem) = MTK.has_ps(getfield(x, :sys))
MTK.get_ps(x::AbstractConductanceSystem) = get_ps(getfield(x, :sys))
MTK.get_defaults(x::AbstractConductanceSystem) = get_defaults(getfield(x, :sys))
MTK.get_states(x::AbstractConductanceSystem) = get_states(getfield(x, :sys))
MTK.get_ivs(x::AbstractConductanceSystem) = get_ivs(getfield(x, :sys))
MTK.get_iv(x::AbstractConductanceSystem) = get_iv(getfield(x, :sys))
MTK.independent_variables(x::AbstractConductanceSystem) = MTK.independent_variables(getfield(x, :sys))
MTK.get_observed(x::AbstractConductanceSystem) = MTK.get_observed(getfield(x, :sys))

function ConductanceSystem(g::Num, ion::IonSpecies, gate_vars::Vector{GatingVariable};
        gbar::Num, linearity::IVCurvature = Linear, extensions::Vector{ODESystem} = ODESystem[],
                           defaults = Dict(), name::Symbol = Base.gensym("Conductance"))

    eqs = Equation[]
    inputs = Set{Num}()
    gate_var_outputs = Set{Num}()
    embed_defaults = Dict()
    params = Set{Num}()
    
    isparameter(gbar) && push!(params, gbar)

    for var in gate_vars
        x, x∞, τₓ = output(var), steadystate(var), timeconstant(var)
        push!(gate_var_outputs, x)
        eq = D(x) ~ inv(τₓ)*(x∞ - x)
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

    return ConductanceSystem(g, ion, gate_vars, inputs, sys, linearity, nothing, name)
end

function IonChannel(ion::IonSpecies,
                    gate_vars::Vector{GatingVariable} = GatingVariable[];
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
                         gate_vars::Vector{GatingVariable} = GatingVariable[],
                         max_s::Union{Num, ElectricalConductance} = 0mS;
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
    gbar_sym = getproperty(cond, :gbar, namespace=false)
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
    gates = GatingVariable[]
    for gate in ex.args
      push!(gates, eval(gate))
    end
    name, ion, gates
end
