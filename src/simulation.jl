"""
$(TYPEDSIGNATURES)

Compile and run a simulation of a single `neuron` or `network` of neurons for a specified
duration, `time`.

If `return_system == true`, returns a simplified `ODESystem` instead.
"""
function Simulation(neuron::AbstractCompartmentSystem; time::Time, return_system = false)
    odesys = convert(ODESystem, neuron)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(odesys)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODAEProblem(simplified, [], (0., t_val), [])
    end
end

function Simulation(neuron::CompartmentSystem{LIF}; time::Time, return_system = false)
    odesys = convert(ODESystem, neuron) & ODESystem([D(@nonamespace(neuron.S_pre)) ~ 0.0],
                                                    t, [], []; name = :foo)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(odesys)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODAEProblem(simplified, [], (0., t_val), [])
    end
end

function Simulation(network::NeuronalNetworkSystem; time::Time, return_system = false)
    odesys = convert(ODESystem, network)
    t_val = ustrip(Float64, ms, time)
    simplified = structural_simplify(odesys)
    if return_system
        return simplified
    else
        @info repr("text/plain", simplified)
        return ODAEProblem(simplified, [], (0., t_val), [])
    end
end


