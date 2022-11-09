
function indexmap(syms, ref)
    idxs = similar(syms, Int64)
    for i in eachindex(syms)
        idxs[i] = indexof(syms[i], ref)
    end
    return idxs
end

simple_spike_check(V,Vprev) = V >= 10. && Vprev < 10.

struct SimpleSpikeDetection
    voltage_indexes::Vector{Int}
end

function SimpleSpikeDetection(network, simplified)
    dvs = states(simplified)
    topo = get_topology(network)
    neuron_Vms = voltage.(neurons(topo))
    Vm_idxs = indexmap(neuron_Vms, dvs)
    return SimpleSpikeDetection(Vm_idxs)
end

function (ssd::SimpleSpikeDetection)(integrator)::BitVector
    idxs = ssd.voltage_indexes
    V = view(integrator.u, idxs)
    Vprev = view(integrator.uprev, idxs)
    return simple_spike_check.(V, Vprev)
end

struct SimpleExponentialResponse
    state_indexes::Vector{Int}
    graph
end

function SimpleExponentialResponse(network, conductance, simplified)
    # calculate state indexes 
    dvs = states(simplified)
    topo = get_topology(network)
    nrns = renamespace.(network, neurons(topo))
    cond_name = nameof(conductance)
    
    for nrn in nrns
        if hasproperty(nrn, cond_name)

        end
    end

    return SimpleExponentialResponse(state_indexes, graph)
end

function simple_spike_propagation!(states, S, graph)
    rows = rowvals(graph) # presynaptic neuron indexes for each stored value
    vals = nonzeros(graph) # weights
    n = size(graph, 1) # n columns
    for i in 1:n
        for j in nzrange(graph, i)
            S[rows[j]] || continue
            states[i] += vals[j]
        end
    end
end

function (ser::SimpleExponentialResponse)(integrator, S)
    idxs = ser.state_indexes
    g = ser.graph
    states = view(integrator.u, idxs)
    spike_propagation_callback!(states, S, g)
    return 
end

struct NetworkCallbacks{D,R,V}
    spike_detection::D
    spike_responses::Vector{R}
    voltage_reset::V
end

function NetworkCallbacks(spike_detection, spike_responses, voltage_reset = identity)
    return NetworkCallbacks(spike_detection, spike_responses, voltage_reset)
end

function (ncs::NetworkCallbacks)(integrator)
    S = ncs.spike_detection(integrator)
    for sr! in spike_responses
        sr!(integrator, S)
    end
    ncs.voltage_reset(integrator)
end

