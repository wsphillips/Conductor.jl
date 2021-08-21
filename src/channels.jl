mutable struct AuxConversion
    params::Vector{Num}
    eqs::Vector{Equation}
end

# Conductance types (conductance as in "g")
abstract type AbstractConductance end

isbuilt(x::AbstractConductance) = x.sys !== nothing

struct IonChannel <: AbstractConductance
    gbar::SpecificConductance # scaling term - maximal conductance per unit area
    conducts::DataType # ion permeability
    inputs::Vector{Num} # cell states dependencies (input args to kinetics); we can infer this
    params::Vector{Num}
    kinetics::Vector{<:AbstractGatingVariable} # gating functions; none = passive channel
    sys::Union{ODESystem, Nothing} # symbolic system
end

function _conductance(::Type{ODESystem},gbar_val::T, gate_vars::Vector{<:AbstractGatingVariable};
                      passive::Bool = false, null_init::Bool = false,
                      name::Symbol) where {T <: Real}

    inputs = Set{Num}()
    states = Set{Num}()
    eqs = Equation[]
    # retrieve all variables present in the RHS of kinetics equations
    # g = total conductance (e.g. g(m,h) ~ ̄gm³h)
    if passive
        params = @parameters g
        defaultmap = Pair[g => gbar_val]
    else
        gates = Set{Num}(getsymbol(x) for x in gate_vars)
        @variables g(t)
        push!(states, g)
        params = @parameters gbar
        defaultmap = Pair[gbar => gbar_val]

        for i in gate_vars
            syms = value.(get_variables(getequation(i)))
            for j in syms
                if j ∉ gates
                    isparameter(j) ? push!(params, j) : push!(inputs, j)
                end
            end
        end
        union!(states, gates, inputs)
        push!(eqs, g ~ gbar * prod(hasexponent(x) ? getsymbol(x)^x.p : getsymbol(x) for x in gate_vars))
        append!(eqs, getequation(x) for x in gate_vars)
        if null_init
            foreach(gate_vars) do x; defaultmap = _merge(defaultmap, Dict(getsymbol(x)=>0.0)) end
        else
            foreach(gate_vars) do x; defaultmap = _merge(defaultmap, getdefaults(x)) end
        end
    end
    system = ODESystem(eqs, t, states, params; defaults = defaultmap, name = name)
    return (collect(inputs), params, system)
end

# General purpose constructor
function IonChannel(conducts::Type{I},
                    gate_vars::Vector{<:AbstractGatingVariable},
                    max_g::SpecificConductance = 0mS/cm^2;
                    passive::Bool = false, name::Symbol) where {I <: Ion}
    # TODO: Generalize to other possible units (e.g. S/F)
    sys_type = passive ? ODESystem : subtype(eltype(gate_vars))
    @assert sys_type <: AbstractSystem "All gate variables must contain an AbstractSystem!"
    @assert ~isequal(sys_type,AbstractSystem) "All gate variables must contain the same AbstracSystem subtype!"

    gbar_val = ustrip(Float64, mS/cm^2, max_g)
    (inputs, params, system) = _conductance(sys_type,gbar_val, gate_vars, passive = passive, name = name)
    return IonChannel(max_g, conducts, inputs, params, gate_vars, system)
end

function (chan::IonChannel)(newgbar::SpecificConductance)
    newchan = deepcopy(chan)
    @set! newchan.gbar = newgbar
    gbar_val = ustrip(Float64, mS/cm^2, newgbar)
    gsym = length(newchan.kinetics) > 0 ?
           value(first(@parameters gbar)) :
           value(first(@parameters g))

    mapping = Dict([gsym => gbar_val])
    new_defaults = merge(defaults(newchan.sys), mapping)
    @set! newchan.sys.defaults = new_defaults
    return newchan
end

# Alias for ion channel with static conductance
function PassiveChannel(conducts::Type{I}, max_g::SpecificConductance = 0mS/cm^2;
                        name::Symbol = Base.gensym(:Leak)) where {I <: Ion}
    gate_vars = AbstractGatingVariable[]
    return IonChannel(conducts, gate_vars, max_g; name = name, passive = true)
end

struct SynapticChannel <: AbstractConductance
    gbar::ElectricalConductance
    conducts::DataType
    reversal::Num
    inputs::Vector{Num}
    params::Vector{Num}
    kinetics::Vector{<:AbstractGatingVariable}
    sys::Union{ODESystem, Nothing}
end

function SynapticChannel(conducts::Type{I}, gate_vars::Vector{<:AbstractGatingVariable},
                         reversal::Num, max_g::ElectricalConductance = 0mS;
                         passive::Bool = false, name::Symbol) where {I <: Ion}

    sys_type = passive ? ODESystem : subtype(eltype(gate_vars))
    @assert sys_type <: AbstractSystem "All gate variables must contain an AbstractSystem!"
    @assert ~isequal(sys_type,AbstractSystem) "All gate variables must contain the same AbstracSystem subtype!"
    gbar_val = ustrip(Float64, mS, max_g)
    (inputs, params, system) = _conductance(sys_type, gbar_val, gate_vars, passive = passive, null_init = true, name = name)
    return SynapticChannel(max_g, conducts, reversal, inputs, params, gate_vars, system)
end

function GapJunction(conducts::Type{I}, reversal::Num, max_g::ElectricalConductance = 0mS;
                     passive::Bool = false, name::Symbol) where {I <: Ion}
    SynapticChannel(conducts, AbstractGatingVariable[], reversal, max_g, passive = true, name)
end

function (chan::SynapticChannel)(newgbar::ElectricalConductance)
    newchan = deepcopy(chan)
    @set! newchan.gbar = newgbar
    gbar_val = ustrip(Float64, mS, newgbar)
    gsym = length(newchan.kinetics) > 0 ?
           value(first(@parameters gbar)) :
           value(first(@parameters g))

    mapping = Dict([gsym => gbar_val])
    new_defaults = merge(defaults(newchan.sys), mapping)
    @set! newchan.sys.defaults = new_defaults
    return newchan
end
