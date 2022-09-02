module Conductor
using DocStringExtensions
using Graphs
using Distributions
using SparseArrays
import SciMLBase
using ModelingToolkit,
      Unitful,
      Unitful.DefaultSymbols,
      InteractiveUtils,
      Symbolics,
      SymbolicUtils,
      Setfield,
      MacroTools

const MTK = ModelingToolkit

import Symbolics:
    get_variables,
    Symbolic,
    value,
    tosymbol,
    VariableDefaultValue,
    wrap, unwrap, Arr,
    scalarize,
    getname

import ModelingToolkit:
    toparam,
    isparameter,
    Equation,
    defaults,
    AbstractSystem,
    get_eqs,
    get_iv,
    get_ivs,
    get_states,
    get_observed,
    get_defaults,
    get_ps,
    get_systems,
    get_continuous_events,
    get_discrete_events,
    get_unit,
    _merge,
    renamespace,
    hasdefault,
    getdefault,
    setdefault,
    AbstractTimeDependentSystem,
    independent_variables,
    get_variables!,
    validate,
    CheckComponents

import ModelingToolkit.SciMLBase: parameterless_type

import Unitful:
    Time,
    TimeUnits,
    Voltage,
    VoltageUnits,
    Current,
    CurrentUnits,
    Molarity,
    ElectricalConductance,
    ElectricalConductanceUnits

import Unitful: mV, mS, cm, µF, mF, µm, pA, nA, mA, µA, ms, mM, µM

import SymbolicUtils:
    FnType,
    symtype,
    operation,
    arguments

import Base: show, display

export Gate, SimpleGate, AlphaBeta, SteadyStateTau, SteadyState, ConstantValue, IonChannel,
       PassiveChannel, SynapticChannel, Synapse, Junction, AxialConductance

export EquilibriumPotential, Equilibrium, Equilibria, MembranePotential, IonCurrent,
       IonConcentration, Concentration, ExtrinsicPotential, Instrinsic, Extrinsic

export MultiCompartmentTopology, NetworkTopology, Population, add_synapse!, add_layer!

export AuxConversion, D
export Simulation
export @named

export Calcium, Sodium, Potassium, Chloride, Cation, Anion, Leak, Ion, NonIonic

export t

export CompartmentSystem, Compartment, ConductanceSystem, Conductance,
       NeuronalNetworkSystem, NeuronalNetwork, MultiCompartmentSystem, MultiCompartment

export output, get_output, timeconstant, steadystate, forward_rate, reverse_rate,
       hasexponent, exponent

export Sphere, Cylinder, Point, Unitless, area, radius, height

export HodgkinHuxley

# Metadata IDs
struct ConductorMaxConductance end

abstract type AbstractConductanceSystem <: AbstractTimeDependentSystem end
abstract type AbstractCompartmentSystem <: AbstractTimeDependentSystem end
abstract type AbstractNeuronalNetworkSystem <: AbstractTimeDependentSystem end

include("util.jl")
include("primitives.jl")
include("gate.jl")
include("conductance.jl")
include("geometry.jl")
include("compartment.jl")
include("multicompartment.jl")
include("network.jl")
include("simulation.jl")
include("populations.jl")
end # module
