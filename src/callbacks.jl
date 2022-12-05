
iseventbased(x) = false
iseventbased(x::ConductanceSystem{<:EventBasedSynapse}) = true
iseventbased(x::CurrentSystem{<:EventBasedSynapse}) = true

discrete_spike_check(V,Vprev)::Bool = V >= 10. && Vprev < 10.

function map_voltage_indices(network, simplified; roots_only = false)
    dvs = states(simplified)
    if roots_only
        comps = root_compartments(get_topology(network))
    else
        comps = compartments(get_topology(network))
    end
    comp_Vms = voltage.(comps)
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
    for affect! in net.spike_affects
        affect!(integrator, i)
    end
end

# with voltage reset
function (net::NetworkAffects{A,T})(integrator, i) where {A,T}
    for affect! in net.spike_affects
        affect!(integrator, i)
    end
    net.tailcall(integrator, i)
end

const SynapticSystem{M} = Union{ConductanceSystem{M},CurrentSystem{M}} where {M<:SynapticModel}

struct SpikeAffect{M}
    synaptic_system::SynapticSystem{M}
    state_indices::Vector{Int}
end

get_system(sa::SpikeAffect) = sa.synaptic_system

############################################################################################

struct ConstantValueEvent{T} <: SummedEventsSynapse
    alpha::T # amount (real value) to perturb by
    state::Num # symbolic state in (each) postsynaptic compartment to perturb
    saturation::T # events ignored when state ≥ saturation
    threshold::T # threshold for AP event
    # weighted_events::Bool # apply weights to alpha
end

function ConstantValueEvent(state::Num; alpha::T = 1.0, saturation::T = floatmax(T), threshold = 10.0mV) where {T}
    thold = ustrip(T, mV, threshold)
    ConstantValueEvent{T}(alpha, state, saturation, thold)
end

function SpikeAffect(synsys::SynapticSystem{T}, network, simplified) where {T<:ConstantValueEvent}
    model = get_model(synsys)
    comps = compartments(get_topology(network))
    state = renamespace(synsys, model.state)
    syn_states = renamespace.(comps, state)
    state_indices = indexmap(syn_states, states(simplified))
    return SpikeAffect{T}(synsys, state_indices)
end

function (cv::SpikeAffect{<:ConstantValueEvent})(integrator, i)
    sys = get_system(cv)
    model = get_model(sys)
    α, sat = model.alpha, model.saturation
    # retrieve model-specific graph layer from multigraph
    g = get_weights(integrator, sys)
    # row numbers of non-zero values from the ith column of the sparse matrix
    rng = nzrange(g,i)
    iszero(length(rng)) && return nothing
    rows = view(rowvals(g), rng) # i.e. indexes of post_synaptic compartment(s)
    for i in view(cv.state_indices, rows)
        integrator.u[i] += ifelse(integrator.u[i] < sat, α, zero(α))
    end
    return nothing
end

