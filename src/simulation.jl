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

spike_check(V,Vprev) = V >= 10 && Vprev < 10

function discrete_spike_detection_callback(integrator, Vm_idxs)::BitArray{1}
    Vm = view(integrator.u, Vm_idxs)
    Vprev = view(integrator.uprev, Vm_idxs)
    return spike_check.(Vm, Vprev)
end

# S is layer independent; graph and Isyn depend on synapse model
function spike_propagation_callback!(Isyn, S, graph)
    rows = rowvals(graph) # presynaptic neuron indexes for each stored value
    vals = nonzeros(graph) # weights
    n = size(graph, 1) # n columns
    for i in 1:n
        for j in nzrange(graph, i)
            S[rows[j]] || continue
            Isyn[i] += vals[j]
        end
    end
end

function generate_network_callbacks()
    dvs = states(simplified)
    ps = parameters(simplified)

    # extract the indexes of somatic voltage
    topo = get_topology(network)
    multigraph = graph(topo)
    neuron_Vms = voltage.(neurons(topo))
    Vm_idxs = indexmap(neuron_Vms, dvs)
    detect_spikes = Base.Fix2(discrete_spike_detection_callback, Vm_idxs)

end

function Simulation(network::NeuronalNetworkSystem; time::Time, return_system = false,
                    jac = false, sparse = false, parallel = nothing)
    t_val, simplified = simplify_simulation(network, time)

    cb = generate_network_callbacks(network, simplified)
    



    # Vm_view = view(u, Vm_idxs) # add to head of callback functions

    #if return_system
    #    return simplified
    #else
    #    @info repr("text/plain", simplified)
    #    return ODEProblem(simplified, [], (0.0, t_val), []; jac, sparse, parallel)
    #end
end
