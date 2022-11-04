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

function Simulation(neuron::CompartmentSystem{LIF}; time::Time, return_system = false,
                    jac = false, sparse = false, parallel = nothing)
    odesys = convert(ODESystem, neuron; with_cb = true)
    t_val, simplified = simplify_simulation(odesys, time)
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
function indexmap(syms, ref)
    idxs = similar(syms, Int64)
    for i in eachindex(syms)
        idxs[i] = indexof(syms[i], ref)
    end
    return idxs
end

spike_check(V,Vprev)::Float64 = V >= 10 && Vprev < 10

function spike_check_callback(integrator, Vm_idxs)
    Vm = view(integrator.u, Vm_idxs)
    Vprev = view(integrator.uprev, Vm_idxs)
    return spike_check.(Vm, Vprev)
end

function spike_propagation_callback(integrator, S)
    topology = topology(integrator.p) # multilayered graph of weights (one layer per synapse type)

    # each synaptic conductance has a parameter value OR cicular buffer of spike input
    # by default the value is zero (no spikes) so we may need a way to reset on the time step after a spike
    view_of_input_ps = view(integrator.p, input_ps_idxs)

    Vth = 10.0 # eventually pull stored parameters for spike threshold

    for layer in keys(mg)
        adjm = mg[layer]
        
        rows = rowvals(adjm) # presynaptic neuron indexes for each stored value
        vals = nonzeros(adjm) # weights
        n = size(adjm, 1) # n columns
        
        for j in 1:n
            for i in nzrange(adjm, j)
                w = vals[i]
                spiked = S[rows[i]] 
                weighted_sum = sum(w[spiked])
                # synaptic_input_parameter_of_neuron_j = weighted_sum
            end
        end
    end
end

function Simulation(network::NeuronalNetworkSystem; time::Time, return_system = false,
                    jac = false, sparse = false, parallel = nothing)
    t_val, simplified = simplify_simulation(network, time)

    # map neuron somatic voltages

    dvs = states(simplified)
    ps = parameters(simplified)
    
    # extract the indexes of somatic voltage
    neuron_Vms = voltage.(neurons(get_topology(network)))
    Vm_idxs = indexmap(neuron_Vms, dvs)
    
    # create an affect!(integrator) that returns a bit vector of spike events
    sc_cb = Base.Fix2(spike_check_callback, Vm_idxs) # function of integrator


    # Vm_view = view(u, Vm_idxs) # add to head of callback functions

    #if return_system
    #    return simplified
    #else
    #    @info repr("text/plain", simplified)
    #    return ODEProblem(simplified, [], (0.0, t_val), []; jac, sparse, parallel)
    #end
end
