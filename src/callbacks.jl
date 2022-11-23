
iseventbased(x) = false
iseventbased(x::ConductanceSystem{<:EventBasedSynapse}) = true
iseventbased(x::CurrentSystem{<:EventBasedSynapse}) = true

function indexmap(syms, ref)
    idxs = similar(syms, Int64)
    for i in eachindex(syms)
        idxs[i] = indexof(syms[i], ref)
    end
    return idxs
end

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

struct SpikeAffect{M}
    model_system::M
    state_indices::Vector{Int}
end

############################################################################################

struct ConstantUpdate{T} <: SummedEventSynapse
    alpha::T
end

function (cu::SpikeAffect{ConstantUpdate})(integrator, i)
    index = cu.state_indices[i]
    alpha = cu.model_system.alpha
    integrator.u[index] += alpha
end

