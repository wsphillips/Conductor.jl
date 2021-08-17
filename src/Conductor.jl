module Conductor

using ModelingToolkit, Unitful, Unitful.DefaultSymbols, InteractiveUtils
using IfElse, Symbolics, SymbolicUtils, Setfield

import Symbolics: get_variables, Symbolic, value, tosymbol, VariableDefaultValue, wrap
import ModelingToolkit: toparam, isparameter, Equation
import SymbolicUtils: FnType

import Unitful: Time, Voltage, Current, Molarity, ElectricalConductance
import Unitful: mV, mS, cm, µF, mF, µm, pA, nA, mA, µA, ms, mM, µM

import Base: show, display

export Gate, AlphaBetaRates, SteadyStateTau, IonChannel, PassiveChannel, SynapticChannel
export EquilibriumPotential, Equilibrium, Equilibria, MembranePotential, MembraneCurrent
export AuxConversion, D, Network
export Soma, Simulation, Concentration, IonConcentration
export @named
export Calcium, Sodium, Potassium, Chloride, Cation, Anion, Leak, Ion

const ℱ = Unitful.q*Unitful.Na # Faraday's constant
const t = let name = :t; only(@parameters $name) end
const D = Differential(t)

# Helper utils -- FIXME: MTK has these now
hasdefault(x::Symbolic) = hasmetadata(x, VariableDefaultValue) ? true : false
hasdefault(x::Num) = hasdefault(ModelingToolkit.value(x))
hasdefault(x) = false    

getdefault(x::Symbolic) = hasdefault(x) ? getmetadata(x, VariableDefaultValue) : nothing
getdefault(x::Num) = getdefault(ModelingToolkit.value(x))

# Basic symbols
function MembranePotential() 
    name = :Vₘ
    return only(@variables $name(t))
end

@enum Location Outside Inside

# Custom Unitful.jl quantities
@derived_dimension SpecificConductance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^3 # conductance per unit area
@derived_dimension SpecificCapacitance 𝐈^2*𝐋^-4*𝐌^-1*𝐓^4 # capacitance per unit area
@derived_dimension ConductancePerFarad 𝐓^-1 # S/F cancels out to 1/s; perhaps find a better abstract type?

# Metadata IDs
abstract type ConductorCurrentCtx end
abstract type ConductorEquilibriumCtx end
abstract type ConductorConcentrationCtx end
abstract type ConductorAggregatorCtx end

# Ion species
abstract type Ion end
abstract type Cation <: Ion end
abstract type Anion <: Ion end
struct Calcium  <: Ion end
struct Sodium <: Ion end
struct Potassium <: Ion end
struct Chloride <: Ion end

const Ca = Calcium
const Na = Sodium
const K = Potassium
const Cl = Chloride
const Mixed = Ion # non-specific ion
const Leak = Mixed

const PERIODIC_SYMBOL = IdDict(Na => :Na, K  => :K, Cl => :Cl, Ca => :Ca, Leak => :l)

# Concentrations of ions
abstract type AbstractConcentration end

struct IonConcentration{I<:Ion, L<:Location, V<:Union{Nothing, Num,Symbolic,Molarity}} <:AbstractConcentration
    ion::Type{I}
    val::V
    loc::L
end

# FIXME: handle default values better
function Concentration(::Type{I}, val = nothing, loc::Location = Inside, name::Symbol = PERIODIC_SYMBOL[I]) where {I <: Ion}
    sym = Symbol(name,(loc == Inside ? "ᵢ" : "ₒ"))
    # FIXME: Not necessarily a parameter when used as a primitive...but we should support
    # this. Use a flag?
    var = #=val isa Molarity ? only(@parameters $sym) :=# only(@variables $sym(t))
    return setmetadata(var,  ConductorConcentrationCtx, IonConcentration(I, val, loc))
end

isconcentration(x::Symbolic) = hasmetadata(x, ConductorConcentrationCtx)
isconcentration(x::Num) = isconcentration(value(x))
getconcentration(x::Symbolic) = isconcentration(x) ? getmetadata(x, ConductorConcentrationCtx) : nothing
getconcentration(x::Num) = getconcentration(value(x))

# Currents
struct MembraneCurrent{I<:Ion,V<:Union{Nothing,Num,Symbolic,Current}}
    ion::Type{I}
    val::V
end

# TODO: add aggregator as field of Membrane current struct
function MembraneCurrent{I}(val = nothing; name::Symbol = PERIODIC_SYMBOL[I], aggregate::Bool = false) where {I <: Ion}
    sym = Symbol("I", name)
    var = val isa Current ? only(@parameters $sym) : only(@variables $sym(t))
    var = setmetadata(var, ConductorCurrentCtx, MembraneCurrent(I, val))
    return setmetadata(var, ConductorAggregatorCtx, aggregate)
end

ismembranecurrent(x::Symbolic) = hasmetadata(x, ConductorCurrentCtx)
ismembranecurrent(x::Num) = ismembranecurrent(ModelingToolkit.value(x))
getmembranecurrent(x::Union{Num, Symbolic}) = ismembranecurrent(x) ? getmetadata(x, ConductorCurrentCtx) : nothing
iontype(x::Union{Num, Symbolic}) = getmembranecurrent(x).ion
isaggregator(x::Union{Num, Symbolic})  = getmetadata(x, ConductorAggregatorCtx)

# Equilibrium potential implicitly defines an ionic gradient
abstract type AbstractIonGradient end

struct EquilibriumPotential{I<:Ion,V<:Union{Num,Symbolic,Voltage}} <: AbstractIonGradient
    ion::Type{I}
    val::V
end

const Equilibrium{I} = EquilibriumPotential{I} 

function EquilibriumPotential{I}(val, name::Symbol = PERIODIC_SYMBOL[I]) where {I <: Ion}
    sym = Symbol("E", name)
    var = val isa Voltage ? only(@parameters $sym) : only(@variables $sym(t))
    return setmetadata(var, ConductorEquilibriumCtx, EquilibriumPotential(I, val))
end

# Alternate constructor
function Equilibria(equil::Vector)
    out = Num[]
    for x in equil
        !(x.first <: Ion) && throw("Equilibrium potential must be associated with an ion type.")
        if typeof(x.second) <: Tuple
            tup = x.second
            typeof(tup[2]) !== Symbol && throw("Second tuple argument for $(x.first) must be a symbol.")
            push!(out, Equilibrium{x.first}(tup...))
        else
            push!(out, Equilibrium{x.first}(x.second))
        end
    end
    return out
end

# Gating variables (as an interface)
abstract type AbstractGatingVariable end

hassteadystate(x::AbstractGatingVariable) = hasfield(typeof(x), :ss) ? !(isnothing(x.ss)) : false
hasexponent(x::AbstractGatingVariable) = hasfield(typeof(x), :p) ? x.p !== one(Float64) : false
getsymbol(x::AbstractGatingVariable) = x.sym
getequation(x::AbstractGatingVariable) = x.df

struct Gate <: AbstractGatingVariable
    sym::Num # symbol/name (e.g. m, h)
    df::Equation # differential equation
    ss::Union{Nothing, Num} # optional steady-state expression for initialization
    p::Float64 # optional exponent (defaults to 1)
end

abstract type AbstractGateModel end
struct SteadyStateTau <: AbstractGateModel end
struct AlphaBetaRates <: AbstractGateModel end

function Gate(::Type{SteadyStateTau}, name::Symbol, ss::Num, tau::Num, p::Real)
    sym = only(@variables $name(t))
    df = D(sym) ~ (ss-sym)/tau # (m∞ - m)/τₘ
    return Gate(sym, df, ss, p)
end

function Gate(::Type{AlphaBetaRates}, name::Symbol, alpha::Num, beta::Num, p::Float64)
    sym = only(@variables $name(t))
    df = D(sym) ~ alpha * (1 - sym) - beta*sym # αₘ(1 - m) - βₘ*m
    ss = alpha/(alpha + beta) # αₘ/(αₘ + βₘ)
    return Gate(sym, df, ss, p)
end

Gate(t::Type{AlphaBetaRates}, name::Symbol, alpha::Num, beta::Num, p::Real) =
Gate(AlphaBetaRates, name, alpha, beta, Float64(p))

# TODO: find a nicer way to do this
function Gate(::Type{SteadyStateTau}; p = one(Float64), kwargs...)
    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:m∞, :τₘ], syms)
        Gate(SteadyStateTau, :m, kwargs[:m∞], kwargs[:τₘ], p)
    elseif issetequal([:h∞, :τₕ], syms)
        Gate(SteadyStateTau, :h, kwargs[:h∞], kwargs[:τₕ], p)
    else
        throw("invalid keywords")
    end
end

function Gate(::Type{AlphaBetaRates}; p = one(Float64), kwargs...)
    syms = keys(kwargs)
    if length(syms) !== 2
        throw("Invalid number of input equations.")
    elseif issetequal([:αₘ, :βₘ], syms)
        Gate(AlphaBetaRates, :m, kwargs[:αₘ], kwargs[:βₘ], p)
    elseif issetequal([:αₕ, :βₕ], syms)
        Gate(AlphaBetaRates, :h, kwargs[:αₕ], kwargs[:βₕ], p)
    elseif issetequal([:αₙ, :βₙ], syms)
        Gate(AlphaBetaRates, :n, kwargs[:αₙ], kwargs[:βₙ], p)
    else
        throw("invalid keywords")
    end
end

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

# Return ODESystem pretty printing for our wrapper types
Base.show(io::IO, ::MIME"text/plain", x::IonChannel) = Base.display(isbuilt(x) ? x.sys : x)

function _conductance(gbar_val::T, gate_vars::Vector{<:AbstractGatingVariable}; 
                      passive::Bool = false, null_init::Bool = false, name::Symbol) where {T <: Real}
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
        append!(defaultmap, getsymbol(x) => (hassteadystate(x) && !null_init) ? x.ss : 0.0 for x in gate_vars) # fallback to zero
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
    gbar_val = ustrip(Float64, mS/cm^2, max_g)
    (inputs, params, system) = _conductance(gbar_val, gate_vars, passive = passive, name = name)
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
    gbar_val = ustrip(Float64, mS, max_g)
    (inputs, params, system) = _conductance(gbar_val, gate_vars, passive = passive, null_init = true, name = name)
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

abstract type Geometry end
abstract type Sphere <: Geometry end
abstract type Cylinder <: Geometry end

struct Compartment{G}
    cap::SpecificCapacitance
    chans::Vector{<:AbstractConductance}
    states::Vector
    params::Vector
    sys::ODESystem
end

function Compartment{Sphere}(channels::Vector{<:AbstractConductance},
                             gradients; name::Symbol, area::Float64 = 0.628e-3, #radius = 20µm,
                             capacitance::SpecificCapacitance = 1µF/cm^2,
                             V0::Voltage = -65mV,
                             holding::Current = 0nA,
                             stimulus::Union{Function,Nothing} = nothing,
                             aux::Union{Nothing, Vector{AuxConversion}} = nothing)
   
    Vₘ = MembranePotential()
    @variables Iapp(t) Isyn(t)
    params = @parameters cₘ aₘ
    grad_meta = getmetadata.(gradients, ConductorEquilibriumCtx) 
    #r_val = ustrip(Float64, cm, radius) # FIXME: make it so we calculate area from dims as needed
 
    systems = []
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
            for x in i.params
                hasdefault(x) && push!(defaultmap, x => getdefault(x))
            end

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
        push!(systems, sys) 
        
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
            if ismembranecurrent(s) && isaggregator(s)
                push!(newstateeqs, s ~ sum(filter(x -> iontype(x) == iontype(s), currents)))
                push!(states, s)
                
            end
        end
        append!(eqs, newstateeqs)
    end
    
    # propagate default parameter values to channel systems
    vm_eq = D(Vₘ) ~ (Iapp - (+(currents..., Isyn)))/(aₘ*cₘ)
    push!(eqs, vm_eq)
    system = ODESystem(eqs, t, states, params; systems = systems, defaults = defaultmap, name = name)
    
    #return (eqs, states, params)
    return Compartment{Sphere}(capacitance, channels, states, params, system)
end

const Soma = Compartment{Sphere}

function Compartment{Cylinder}() end
# should also be able to parse "collections" of compartments" that have an adjacency list/matrix 

# takes a topology, which for now is just an adjacency list; also list of neurons, but we
# should be able to just auto-detect all the neurons in the topology
function Network(neurons, topology; name = :Network)
    
    all_neurons = Set(getproperty.(neurons, :sys))
    eqs = Equation[]
    params = Set()
    states = Set() # place holder
    defaultmap = Pair[]
    systems = []
    push!(systems, all_neurons...)
    post_neurons = Set()

    # what types of synapses do we have
    all_synapses = [synapse.second[2] for synapse in topology]
    all_synapse_types = [x.sys.name for x in all_synapses]
    synapse_types = unique(all_synapse_types)

    # how many of each synapse type are there 
    synapse_counts = Dict{Symbol,Int64}()

    for n in synapse_types
        c = count(x -> isequal(x, n), all_synapse_types)
        push!(synapse_counts, n => c)
    end
    
    # Extract reversal potentials
    reversals = unique([x.reversal for x in all_synapses])
    push!(params, reversals...)
    rev_meta = [getmetadata(x, ConductorEquilibriumCtx).val for x in reversals]
    for (i,j) in zip(reversals, rev_meta)
    push!(defaultmap, i => ustrip(Float64, mV, j))
    end
    
    voltage_fwds = Set()

    # create a unique synaptic current for each post-synaptic target
    for synapse in topology
        pre = synapse.first.sys # pre-synaptic neuron system
        post = synapse.second[1].sys # post-synaptic neuron system
        push!(post_neurons, post)
        Erev = synapse.second[2].reversal
        syntype = synapse.second[2].sys # synapse system
        syn = @set syntype.name = Symbol(syntype.name, synapse_counts[syntype.name]) # each synapse is uniquely named
        synapse_counts[syntype.name] -= 1 
        push!(systems, syn)
        push!(voltage_fwds, syn.Vₘ ~ pre.Vₘ)
        
        if post.Isyn ∈ Set(vcat((get_variables(x.lhs) for x in eqs)...))
            idx = findfirst(x -> isequal([post.Isyn], get_variables(x.lhs)), eqs) 
            eq = popat!(eqs, idx)
            eq = eq.lhs ~ eq.rhs + (syn.g * (post.Vₘ - Erev))
            push!(eqs, eq)
        else
            push!(eqs, post.Isyn ~ syn.g * (post.Vₘ - Erev))
        end
    end

    for nonpost in setdiff(all_neurons, post_neurons)
        push!(eqs, D(nonpost.Isyn) ~ 0)
    end
    
    append!(eqs, collect(voltage_fwds))
    return ODESystem(eqs, t, states, params; systems = systems, defaults = defaultmap, name = name )
end

function Simulation(network; time::Time)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(network)
    display(simplified)
    return ODAEProblem(simplified, [], (0., t_val), [])
end

function Simulation(neuron::Soma; time::Time)
    system = neuron.sys
    # for a single neuron, we just need a thin layer to set synaptic current constant
    @named simulation = ODESystem([D(system.Isyn) ~ 0]; systems = [system])
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(simulation)
    display(simplified)
    return ODAEProblem(simplified, [], (0., t_val), [])
end

end # module

