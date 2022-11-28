
iseventbased(x) = false
iseventbased(x::ConductanceSystem{<:EventBasedSynapse}) = true
iseventbased(x::CurrentSystem{<:EventBasedSynapse}) = true

discrete_spike_check(V,Vprev)::Bool = V >= 10. && Vprev < 10.

function map_voltage_indices(network, simplified)
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
    voltage_indices::Vector{Int}
end

# checks all neurons for threshhold crossing, storing the result in `out`
function (csd::ContinuousSpikeDetection)(out, u, t, integrator)
    for (i, j) in enumerate(csd.voltage_indices)
        out[i] = u[j] - 10.0
    end
    return nothing
end

# voltage reset optional; determined by _compartment model_ not synapse
struct VoltageReset{T}
    voltage_indices::Vector{Int}
    V_rest::T
end

function (vr::VoltageReset)(integrator, i)
    idx = vr.voltage_indices[i]
    integrator.u[idx] = vr.V_rest
end

struct NetworkAffects{A,T}
    spike_affects::Vector{A}
    tailcall::T
end

function NetworkAffects(spike_affects, tailcall = nothing)
    return NetworkAffects(spike_affects, tailcall)
end

# default when we have no voltage reset
function (net::NetworkAffects{A,Nothing})(integrator, i) where {A}
    for affect! in ncs.spike_affects
        affect!(integrator, i)
    end
end

# with voltage reset
function (net::NetworkAffects{A,T})(integrator, i) where {A,T}
    for affect! in ncs.spike_affects
        affect!(integrator, i)
    end
    ncs.tailcall(integrator, i)
end

const SynapticSystem{M} = Union{ConductanceSystem{M},CurrentSystem{M}} where {M<:SynapticModel}

struct SpikeAffect{M}
    synaptic_system::SynapticSystem{M}
    state_indices::Vector{Int}
end

get_model(sa::SpikeAffect) = get_model(sa.synaptic_system)

############################################################################################

struct ConstantValueEvent{T} <: SummedEventSynapse
    alpha::T # amount (real value) to perturb by
    state::Num # symbolic state in (each) postsynaptic compartment to perturb
    # weighted_events::Bool # apply weights to alpha or not
    # saturation::T # events ignored when state ≥ saturation
end

function SpikeAffect(synsys::SynapticSystem{ConstantValueEvent}, network, simplified)
    model = get_model(synsys)
    comps = compartments(topology(network))
    syn_states = renamespace.(comps, model.state)
    state_indices = indexmap(syn_states, states(simplified))
    return SpikeAffect{ConstantValueEvent}(synsys, state_indices)
end

function (cv::SpikeAffect{ConstantValueEvent})(integrator, i)
    model = get_model(cv)
    # retrieve model-specific graph layer from multigraph
    g = get_weights(integrator, model)
    # row numbers of non-zero values from the ith column of the sparse matrix
    rows = view(rowvals(g), nzrange(g,i)) # i.e. indexes of post_synaptic compartment(s)
    for i in view(cv.state_indices, rows)
        integrator.u[i] += model.alpha
    end
end

struct WeightedValueEvent{T} <: SummedEventSynapse
    alpha::T
    state::Num
end
