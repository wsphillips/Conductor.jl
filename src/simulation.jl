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
                    jac = false, sparse = false, parallel = nothing)
    t_val, simplified = simplify_simulation(network, time)
    if return_system
        return simplified
    else
        if any(x -> typeof(x) <: EventBasedSynapse, keys(graph(topology(network))))
            cb_affect = generate_network_callback(network, simplified)
            dcb = DiscreteCallback(CONSTANT_CONDITION, cb_affect)
            ODEProblem(simplified, [], (0.0, t_val), [];
                       callback = dcb, jac, sparse, parallel)
        else
            ODEProblem(simplified, [], (0.0, t_val), []; jac, sparse, parallel)
        end
    end
end

function generate_network_callback(network, simplified)
    spike_detection =
    spike_affects = []

    for _ in XXX
        push!(spike_affects) = YYY  
    end

    tailcall = identity
    return NetworkCallbacks(spike_detection, spike_affects, tailcall)
end
