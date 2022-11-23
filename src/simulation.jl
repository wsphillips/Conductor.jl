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
function Simulation(neuron::AbstractCompartmentSystem; time::Time, return_system = false,
                    jac = false, sparse = false,
                    parallel = nothing)
    t_val, simplified = simplify_simulation(neuron, time)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODEProblem(simplified, [], (0.0, t_val), []; jac, sparse, parallel)
    end
end

struct NetworkParameters{T,G,C}
    ps::Vector{T}
    topology::NetworkTopology{G,C}
end

Base.getindex(x::NetworkParameters, i) = x.ps[i]
topology(x::NetworkParameters) = getfield(x, :topology)
indexof(sym, syms) = findfirst(isequal(sym), syms)

function Simulation(network::NeuronalNetworkSystem; time::Time, return_system = false,
                    jac = false, sparse = false, parallel = nothing, continuous_events = false)
    t_val, simplified = simplify_simulation(network, time)
    return_system && return simplified
    event_models = filter(x -> typeof(x) <: EventBasedSynapse, synaptic_models(network))
    if isempty(event_models) # ie. all continuous synapses
        return ODEProblem(simplified, [], (0.0, t_val), []; jac, sparse, parallel)
    else
        for model in event_models
            generate_callback(model, network, simplified; continuous_events)
        end
        ODEProblem(simplified, [], (0.0, t_val), []; callback = cb, jac, sparse, parallel)
    end
end

# if continuous, we use a vector condition: cond(out, u, t, integrator)
# conditions should be synaptic model agnostic, but customizable (defines what a spike is)
function generate_callback_condition(network, simplified; continuous_events)
    Vm_idxs = map_voltage_indices(network, simplified)
    if continuous_events
        # returns a single functor of form cond(out, u, t, integrator)::Vector{eltype(u)}
        return ContinuousSpikeDetection(Vm_idxs)
    else # n (# neurons) discrete conditions: cond(u, t, integrator)::Bool 
        return [DiscreteSpikeDetection(x) for x in Vm_idxs]
    end
end

function generate_callback_affect(network, simplified; continuous_events)
    spike_affects = Dict() 
    if continuous_events
        # return single functor affect compatible with `VectorContinuousCallback`
        # `affect(integrator, i)` where i denotes a spike from the i-th neuron
        # applies the spike response for each synaptic model/network layer sequentially
        for model in layers(network)
            spike_affects[model] = SpikeAffect(model, network, simplified)
        end
    else
        # return vector of affects
        # return two argument functor `affect(integrator, i)` as above, but we will apply
        # `Base.Fix2(integrator, i)` for each neuron a priori. In the discrete case, spike
        # event conditions and affects are applied independently in a long `CallbackChain`

        nrn_affects = []
        for model in layers(network)
            for (i, neuron) in neurons(network)
                push!(nrn_affects,
                      Base.Fix2(SpikeAffect(model, network, simplified), i))
            end
        end
        spike_affects[model] = nrn_affects
    end
    tailcall = identity # placeholder for voltage reset
    return NetworkCallbacks(spike_detection, spike_affects, tailcall)
end

function generate_callback(model, network, simplified; continuous_events)
    cb_condition = generate_callback_condition(network, simplified; continuous_events)
    cb_affect = generate_callback_affect(network, simplified; continuous_events)
    if continuous_events
        return VectorContinuousCallback(cb_condition, cb_affect) 
    else
        callbacks = DiscreteCallback.(cb_condition, cb_affect)   
        return CallbackChain(callbacks...)
    end
end
