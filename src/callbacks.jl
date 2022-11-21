
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

struct SymbolicSpikeAffect{F,S}
    state_indexes::Vector{Int} # the index of the modified state in every compartment
    update_indexes::Vector{Tuple{Vector{Int},Vector{Int}}} # stencils for dvs and ps used by update function
    model_system::S
    update::F # function that computes the amount to perturb the synaptic state
end

function SymbolicSpikeAffect(model_system::SymbolicUpdateSynapse, network, simplified)
    # calculate state indexes 
    ###########################################
    dvs = states(simplified) # symbols in order of lowered system
    ps = parameters(simplified)
    rule = update_rule(model_system)
    updated_state = rule.lhs
    update_expr = rule.rhs

    # symbol for state namespaced wrt the conductance model
    tmp = renamespace(model_system, updated_state)
    
    # fetches the state from every compartment
    states_to_update = []
    for sys in systems(network)
        push!(states_to_update, renamespace(sys, tmp))
    end

    state_indexes = indexmap(states_to_update, dvs)
    ############################################
    # update rule RGF

    update_dvs, update_ps = [], []
    collect_vars!(update_dvs, update_ps, update_expr, t)
    
    # returns RGF of form: update(dvs, ps) -> (new value of synaptic state)
    update_function = @RuntimeGeneratedFunction(build_function(update_expr, update_dvs, update_ps))
    
    # renamespace to the current/conductance system
    update_dvs, update_ps = renamespace.(model_system, update_dvs), renamespace.(model_system, update_ps)

    # vector of tuples containing the index stencil for dvs and ps of each compartment
    update_idx_pairs = []
    for sys in systems(network)
        idxs_dvs = indexmap(renamespace.(sys, update_dvs), dvs)
        idxs_ps = indexmap(renamespace.(sys, update_ps), ps)
        push!(update_idx_pairs, (idxs_dvs, idxs_ps))            
    end

    return SpikeAffect(state_indexes, update_idx_pairs, model_system, update_function)
end

# for discrete callbacks, use Base.Fix2(sa, i)
# weights matrix should be a column-wise presynaptic neurons (transpose of standard form)
# this is a sparse update; it might be simpler (but potentially expensive) to do dense?
function (sa::SymbolicSpikeAffect{F,S})(integrator, i)
    g = weights(integrator.p)[sa.model_system] # layer-specific adjacency matrix
    update = sa.update 
    rows = rowvals(g)
    vals = nonzeros(g)
    for j in nzrange(g, i) # index range of the nonzero values in the ith column
        post = rows[j] # each post synaptic target
        W = vals[j] # weight of connection
        state_index, update_indexes = sa.state_indexes[post], sa.update_indexes[post]
        dvs, ps = view(integrator.u, update_indexes[1]), view(integrator.p, update_indexes[2])
        integrator.u[post] += update(dvs, ps) * W
    end
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
#=
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
=#

