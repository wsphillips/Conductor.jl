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

CONSTANT_CONDITION(u, t, integrator) = true

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
# conditions are model agnostic
function generate_callback_condition(network, simplified; continuous_events)
    if continuous_events
       # return functor condition compatible with VectorContinuousCallback 
       # cond(out, u, t, integrator)::Vector{eltype(u)}
       # where out[i] == 0.0 for each event
       # each neuron -> each event
    else
       # n (# neurons) discrete conditions: cond(u, t, integrator)::Bool 
    end
end

function generate_callback_affect(network, simplified; continuous_events)
    if continuous_events
        # return functor affect compatible with VectorContinuousCallback
        # affect(integrator, i) where i denotes a spike from the i-th neuron
        # applies the spike response for each synaptic model/network layer sequentially
    else
        # return functor  
    end

    tailcall = identity # placeholder for 
    return NetworkCallbacks(spike_detection, spike_affects, tailcall)
end

function generate_callback(model, network, simplified; continuous_events)
    cb_condition = generate_callback_condition(network, simplified; continuous_events)
    cb_affect = generate_callback_affect(network, simplified; continuous_events)
    return CallbackChain()   
end
