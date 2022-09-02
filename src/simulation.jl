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
        return ODEProblem(simplified, [], (0., t_val), []; jac, sparse, parallel)
    end
end

function Simulation(neuron::CompartmentSystem{LIF}; time::Time, return_system = false,
                                                    jac = false, sparse = false, parallel = nothing)
    odesys = convert(ODESystem, neuron; with_cb = true)
    t_val, simplified = simplify_simulation(odesys, time)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODEProblem(simplified, [], (0., t_val), []; jac, sparse, parallel)
    end
end

function Simulation(network::NeuronalNetworkSystem; time::Time, return_system = false,
                                                    jac = false, sparse = false, parallel = nothing)
    t_val, simplified = simplify_simulation(network, time)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODEProblem(simplified, [], (0., t_val), []; jac, sparse, parallel)
    end
end


