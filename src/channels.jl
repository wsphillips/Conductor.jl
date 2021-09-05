# Conductance types (conductance as in "g")
abstract type AbstractConductanceSystem <: AbstractTimeDependentSystem end

# Linear Ohmic/Nernst vs non-linear GHK
@enum IVCurvature Linear Rectifying

# TODO: Implement AbstractSystem methods for AbstractConductanceSystem
getion(x::AbstractConductanceSystem) = getfield(x, :ion)
getinputs(x::AbstractConductanceSystem) = getfield(x, :inputs)
# Abstract types without parametrics
struct ConductanceSystem{S<:AbstractTimeDependentSystem} <: AbstractConductanceSystem
    output::Num # 'g' by default 
    ion::IonSpecies # ion permeability
    gate_vars::Vector{<:AbstractGatingVariable}
    inputs::Vector # required inputs
    sys::S
    linearity::IVCurvature
    transmatrix::Union{Matrix, Nothing}
end

function ConductanceSystem(g::Num, ion::IonSpecies, gate_vars::Vector{GatingVariable};
                           max_g::Real = 0.0, linearity::IVCurvature = Linear,
                           defaults = Dict(), name::Symbol = Base.gensym("Conductance"))

    eqs = Equation[]
    inputs = Set{Num}()
    gate_var_outputs = Set{Num}()
    embed_defaults = Dict()
    params = Set{Num}(@parameters gbar)
    push!(defaults, gbar => max_g)

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

    setdiff!(inputs, params, gate_var_outputs)

    push!(eqs, g ~ gbar * prod(hasexponent(x) ? output(x)^exponent(x) : output(x) for x in gate_vars))
    sys = ODESystem(eqs, t, union(inputs, g, gate_var_outputs), params;
                    defaults = merge(embed_defaults, defaults), name = name)

    return ConductanceSystem(g, dimension(max_g), ion, gate_vars, sys, linearity, transmatrix, nothing)
end

function IonChannel(ion::IonSpecies,
                    gate_vars::Vector{GatingVariable} = GatingVariable[];
                    max_g::SpecificConductance = 0mS/cm^2,
                    name::Symbol = Base.gensym("IonChannel"),
                    linearity::IVCurvature = Linear, defaults = Dict())

    @variables g(t)
    gbar_val = ustrip(Float64, mS/cm^2, max_g)
    setmetadata(g, ConductorUnits, mS/cm^2) # TODO: rework with MTK native unit system
    ConductanceSystem(g, ion, gate_vars;
                      max_g = gbar_val, name = name, defaults = defaults, linearity = linearity)
end

function SynapticChannel(ion::IonSpecies,
                         gate_vars::Vector{GatingVariable} = GatingVariable[],
                         max_s::ElectricalConductance = 0mS;
                         name::Symbol = Base.gensym("SynapticChannel"),
                         linearity::IVCurvature = Linear, defaults = Dict())
    @variables s(t)
    sbar_val = ustrip(Float64, mS, max_s)
    setmetadata(s, ConductorUnits, mS) # TODO: rework iwth MTK native unit system
    ConductanceSystem(s, ion, gate_vars;
                      max_g = sbar_val, name = name, defaults = defaults, linearity = linearity)
end

function (cond::AbstractConductanceSystem)(newgbar::Quantity)
    newcond = deepcopy(cond)
    g = output(newcond)
    outunits = getmetadata(g, ConductorUnits)
    if dimension(outunits) !== dimension(newgbar)
        @error "Input Dimensions do not match output of ConductanceSystem"
    end
    gbar_val = ustrip(Float64, outunits, newgbar)
    gbar_sym = getproperty(cond, gbar, namespace=false)
    push!(defaults(newcond), gbar_sym => gbar_val)
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
