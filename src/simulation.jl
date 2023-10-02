
# for now, fallback to ODESystem conversion
function ModelingToolkit.ODESystem(sys::AbstractConductorSystem; simplify = true)
    odesys = convert(ODESystem, sys)
    return simplify ? structural_simplify(odesys) : odesys
end

# simplified = settunable(simplified, tunable)

"""
$(TYPEDSIGNATURES)

Compile and run a simulation of a single `neuron` or `network` of neurons for a specified
duration, `time`.

If `return_system == true`, returns a simplified `ODESystem` instead.
"""
function Simulation(neuron::AbstractCompartmentSystem, tspan; fixed_params = [],
                    simplify = true, parallel = Symbolics.SerialForm(), kwargs...)
    if isempty(fixed_params)
        simplified = ODESystem(neuron; simplify)
    else
        odesys = convert(ODESystem, neuron)
        simplified = set_fixed_params(odesys, fixed_params; simplify) 
    end
    tstart, tstop = time_span(tspan)
    return ODEProblem(simplified, [], (tstart, tstop), []; parallel, kwargs...)
end

struct NetworkParameters{T}
    ps::Vector{T}
    topology::NetworkTopology
end

Base.getindex(x::NetworkParameters, i) = x.ps[i]
topology(x::NetworkParameters) = getfield(x, :topology)

function get_weights(integrator, model)
    topo = topology(integrator.p)
    return graph(topo)[model]
end

function Simulation(network::NeuronalNetworkSystem, tspan, cb = nothing; simplify = true,
                    parallel = Symbolics.SerialForm(), continuous_events = false,
                    refractory = true, kwargs...)
    odesys = ODESystem(network; simplify)
    tstart, tstop = time_span(tspan)
    if !any(iseventbased.(synaptic_systems(network)))
        return ODEProblem(odesys, [], (tstart, tstop), []; parallel, kwargs...)
    else
        if isnothing(cb)
            cb = generate_callback(network, odesys; continuous_events, refractory)
        end
        prob = ODEProblem(odesys, [], (tstart, tstop), [];
                          callback = cb, parallel, kwargs...)
        remake(prob; p = NetworkParameters(prob.p, get_topology(network)))
    end
end

# if continuous, condition has vector cb signature: cond(out, u, t, integrator)
function generate_callback_condition(network, odesys; continuous_events, refractory)
    voltage_indices = map_voltage_indices(network, odesys; roots_only = true)
    if continuous_events
        return ContinuousSpikeDetection(voltage_indices)
    else # discrete condition for each compartment
        return [DiscreteSpikeDetection(voltage_index, refractory) for voltage_index in voltage_indices]
    end
end

function generate_callback_affects(network, odesys)
    spike_affects = []
    for sys in synaptic_systems(network)
        push!(spike_affects, SpikeAffect(sys, network, odesys))
    end
    tailcall = nothing # placeholder for voltage reset
    return NetworkAffects(spike_affects, tailcall)
end

function generate_callback(network, odesys; continuous_events, refractory)
    cb_condition = generate_callback_condition(network, odesys; continuous_events, refractory)
    cb_affect = generate_callback_affects(network, odesys)
    if continuous_events
        return VectorContinuousCallback(cb_condition, cb_affect,
                                        length(cb_condition.voltage_indices)) 
    else
        affects = []
        for i in 1:length(root_compartments(get_topology(network)))
            push!(affects, Base.Fix2(cb_affect, i))
        end
        callbacks = [DiscreteCallback(x,y) for (x,y) in zip(cb_condition, affects)]
        return CallbackSet(callbacks...)
    end
end
