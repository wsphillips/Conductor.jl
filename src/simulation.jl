function simplify_simulation(sys, time)
    odesys = convert(ODESystem, sys)
    t_val = ustrip(Float64, ms, time)
    return t_val, structural_simplify(odesys)
end

"""
$(TYPEDSIGNATURES)

Compile and run a simulation of a single `neuron` or `network` of neurons for a specified
duration, `time`.

If `return_system == true`, returns a simplified `ODESystem` instead.
"""
function Simulation(neuron::AbstractCompartmentSystem, time::Time; return_system = false,
                    jac = false, sparse = false,
                    parallel = Symbolics.SerialForm())
    t_val, simplified = simplify_simulation(neuron, time)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODEProblem(simplified, [], (0.0, t_val), []; jac, sparse, parallel)
    end
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
function Simulation(network::NeuronalNetworkSystem, time::Time; return_system = false,
        jac = false, sparse = false, parallel = Symbolics.SerialForm(), continuous_events = false)
    t_val, simplified = simplify_simulation(network, time)
    return_system && return simplified
    if !any(iseventbased.(synaptic_systems(network)))
        return ODEProblem(simplified, [], (0.0, t_val), []; jac, sparse, parallel)
    else
        cb = generate_callback(network, simplified; continuous_events)
        prob = ODEProblem(simplified, [], (0.0, t_val), []; callback = cb, jac, sparse, parallel)
        remake(prob; p = NetworkParameters(prob.p, get_topology(network) ))
    end
end

# if continuous, condition has vector cb signature: cond(out, u, t, integrator)
function generate_callback_condition(network, simplified; continuous_events)
    voltage_indices = map_voltage_indices(network, simplified; roots_only = true)
    if continuous_events
        return ContinuousSpikeDetection(voltage_indices)
    else # discrete condition for each compartment
        return [DiscreteSpikeDetection(voltage_index) for voltage_index in voltage_indices]
    end
end

function generate_callback_affects(network, simplified)
    spike_affects = []
    for sys in synaptic_systems(network)
        push!(spike_affects, SpikeAffect(sys, network, simplified))
    end
    tailcall = nothing # placeholder for voltage reset
    return NetworkAffects(spike_affects, tailcall)
end

function generate_callback(network, simplified; continuous_events)
    cb_condition = generate_callback_condition(network, simplified; continuous_events)
    cb_affect = generate_callback_affects(network, simplified)
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
