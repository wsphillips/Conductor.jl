# takes a topology, which for now is just an adjacency list; also list of neurons, but we
# should be able to just auto-detect all the neurons in the topology

abstract type AbstractNetworkSystem <: AbstractTimeDependentSystem

struct NetworkSystem <: AbstractNetworkSystem end
    topology
    extensions
    defaults
    name
end

function NetworkSystem(neurons, topology; name = Base.gensym(:Network))

    all_neurons = Set(getproperty.(neurons, :sys))
    eqs = Equation[]
    params = Set()
    states = Set() # place holder
    defaultmap = Pair[]
    all_systems = AbstractSystem[]
    push!(all_systems, all_neurons...)
    post_neurons = Set()

    # what types of synapses do we have
    all_synapses = [synapse.second[2] for synapse in topology]
    all_synapse_types = [nameof(x.sys) for x in all_synapses]
    synapse_types = unique(all_synapse_types)

    # how many of each synapse type are there
    synapse_counts = Dict{Symbol,Int64}()

    for n in synapse_types
        c = count(x -> isequal(x, n), all_synapse_types)
        push!(synapse_counts, n => c)
    end

    # Extract reversal potentials
    reversals = unique([x.reversal for x in all_synapses])
    push!(params, reversals...)
    rev_meta = [getmetadata(x, ConductorEquilibriumCtx).val for x in reversals]

    for (i,j) in zip(reversals, rev_meta)
        push!(defaultmap, i => ustrip(Float64, mV, j))
    end

    voltage_fwds = Set()

    # create a unique synaptic current for each post-synaptic target
    for synapse in topology
        pre = synapse.first.sys # pre-synaptic neuron system
        post = synapse.second[1].sys # post-synaptic neuron system
        push!(post_neurons, post)
        Erev = synapse.second[2].reversal
        syntype = synapse.second[2].sys # synapse system
        synname = nameof(syntype)
        syn = @set syntype.name = Symbol(synname, synapse_counts[synname]) # each synapse is uniquely named
        synapse_counts[synname] -= 1
        push!(all_systems, syn)
        push!(voltage_fwds, syn.Vₘ ~ pre.Vₘ)

        if post.Isyn ∈ Set(vcat((get_variables(x.lhs) for x in eqs)...))
            idx = findfirst(x -> isequal([post.Isyn], get_variables(x.lhs)), eqs)
            eq = popat!(eqs, idx)
            eq = eq.lhs ~ eq.rhs + (syn.g * (post.Vₘ - Erev))
            push!(eqs, eq)
        else
            push!(eqs, post.Isyn ~ syn.g * (post.Vₘ - Erev))
        end
    end

    for nonpost in setdiff(all_neurons, post_neurons)
        push!(eqs, D(nonpost.Isyn) ~ 0)
    end

    append!(eqs, collect(voltage_fwds))
    network_system = ODESystem(eqs, t, states, params; defaults = defaultmap, name = name )
    return compose(network_system, all_systems)
end
