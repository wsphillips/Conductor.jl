
function indexmap(syms, ref)
    idxs = similar(syms, Int64)
    for i in eachindex(syms)
        idxs[i] = indexof(syms[i], ref)
    end
    return idxs
end

discrete_spike_check(V,Vprev)::Bool = V >= 10. && Vprev < 10.

function voltage_indexes(network, simplified)
    dvs = states(simplified)
    compartments = get_systems(network)
    comp_Vms = renamespace.(network, voltage.(compartments))
    return indexmap(comp_Vms, dvs)
end

struct DiscreteSpikeDetection
    voltage_index::Int
end

# checks whether a pre-determined neuron spiked
function (dsd::DiscreteSpikeDetection)(u, t, integrator)
    idx = dsd.voltage_index
    V = u[idx]
    Vprev = integrator.uprev[idx]
    return discrete_spike_check(V, Vprev)
end

struct ContinuousSpikeDetection
    voltage_indexex::Vector{Int}
end

# checks all neurons for threshhold crossing, storing the result in `out`
function (csd::ContinuousSpikeDetection)(out, u, t, integrator)
    for (i, j) in enumerate(csd.voltage_indexes)
        out[i] = u[j] - 10.0
    end
    return nothing
end

struct SpikeAffect{F,S}
    state_indexes::Vector{Int} # the index of the modified state in every compartment
    quantity_indexes::Vector{Tuple{Vector{Int},Vector{Int}}}
    model_system::S
    quantity::F # function that computes the 
end

function SpikeAffect(model_system, network, simplified)
    # calculate state indexes 
    dvs = states(simplified) # symbols in order of lowered system
    ps = parameters(simplified)
    rule = update_rule(model_system)

    # symbol for state namespaced wrt the conductance model
    tmp = renamespace(model_system, updated_state)
    
    # fetches the state from every compartment
    states_to_update = []
    for sys in systems(network)
        push!(states_to_update, renamespace(network, renamespace(sys, tmp)))
    end

    state_indexes = indexmap(states_to_update, dvs)
    ############################################
    updated_state = rule.lhs
    update_quantity = rule.rhs

    q_dvs, q_ps = [], []
    collect_vars!(q_dvs, q_ps, update_quantity, t)

    calc_quantity = @RuntimeGeneratedFunction(build_function(update_quantity, q_dvs, q_ps))
    
    # we have to map the indexes of the states and parameters used by the quantity calculation
    # and order them the same as our compartments
    # if the model is current based, we scale the synaptic weight by the update quantity

    # 
    q_dvs, q_ps = renamespace.(model_system, q_dvs), renamespace.(model_system, q_ps)
    q_idx_pairs = []
    for sys in systems(network)
        idxs_dvs = indexmap(renamespace.(sys, q_dvs), dvs)
        idxs_ps = indexmap(renamespace.(sys, q_ps), ps)
        push!(q_idx_pairs, (indxs_dvs, idxs_ps))            
    end

    return SpikeAffect(state_indexes, q_idx_pairs, model_system, calc_quantity)
end

# for discrete callbacks, use Base.Fix2(sa, i)
function (sa::SpikeAffect{F,S})(integrator, i) where {S<:ConductanceSystem{SummedEventsSynapse}}
    state_idxs = sa.state_indexes
    q_idx_pairs = quantity_indexes
    multigraph = weights(integrator.p) # needs to be correct layer for model
    g = multigraph[sa.model_system]
    states = view(integrator.u, idxs)
    q_ps, q_dvs =  
    

    return 
end

struct NetworkCallbacks{A,T}
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


