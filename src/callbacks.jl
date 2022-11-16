
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
    comp_Vms = renamespace.(network, voltage.(vertices(topo)))
    Vm_idxs = indexmap(comp_Vms, dvs)
    return SimpleSpikeDetection(Vm_idxs)
end

function (ssd::SimpleSpikeDetection)(integrator)::BitVector
    idxs = ssd.voltage_indexes
    V = view(integrator.u, idxs)
    Vprev = view(integrator.uprev, idxs)
    return simple_spike_check.(V, Vprev)
end

struct SpikeAffect{G,F}
    state_indexes::Vector{Int}
    graph::G
    affect!::F
end

function SpikeAffect(conductance, network, simplified)
    SpikeAffect()
end

function SpikeAffect(conductance, network, simplified)
    # calculate state indexes 
    dvs = states(simplified) # symbols in order of lowered system
    local_states = UpdatedStates(conductance) # local variable names in synapse

    syss = get_systems(network)
    compartments = renamespace.(network, syss)

    cond_name = nameof(conductance)
    
    for comp in compartments
        if hasproperty(comp, cond_name) # sanity check if the compartment has our syn model
        end
    end
    
    state_indexes = indexmap()

    return SpikeAffect(state_indexes, graph, affect!)
end

function (sa::SpikeAffect)(integrator, S)
    idxs = sa.state_indexes
    g = sa.graph
    states = view(integrator.u, idxs)
    sa.affect!(states, S, g)
    return 
end

struct NetworkCallbacks{D,A,T}
    spike_detection::D
    spike_affects::Vector{A}
    tailcall::T
end

function NetworkCallbacks(spike_detection, spike_affects, tailcall = identity)
    return NetworkCallbacks(spike_detection, spike_affects, tailcall)
end

function (ncs::NetworkCallbacks)(integrator)
    S = ncs.spike_detection(integrator)
    for affect! in ncs.spike_affects
        affect!(integrator, S)
    end
    ncs.tailcall(integrator)
end

############################################

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


