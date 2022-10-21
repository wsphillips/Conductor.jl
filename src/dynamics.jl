abstract type AbstractDynamics end

struct HodgkinHuxley <: AbstractDynamics
    "Ionic conductances."
    channels::Vector{AbstractConductanceSystem}
    "Equilibrium potentials belonging to ionic membrane conductances."
    channel_reversals::Vector{Num}
end

function process_reversals!(gen, dynamics::HodgkinHuxley, synapses, arbor)
    (; dvs, ps, eqs) = gen
    reversal_equation_vars = Set{Num}()
    all_reversals = dynamics.channel_reversals # channel reversal potentials
    append!(all_reversals, reversal.(synapses)) # synaptic reversal potentials
    reversals!(all_reversals, arbor) # axial reversal potentials

    for Erev in all_reversals
        if isparameter(Erev)
            push!(ps, Erev)
        else
            push!(dvs, Erev)
            # FIXME: hacked solution. This assumes the equation for a dynamic reversal == default value
            if MTK.hasdefault(Erev)
                get_variables!(reversal_equation_vars, getdefault(Erev))
                push!(eqs, Erev ~ getdefault(Erev))
            end
        end
    end
    filter!(x -> !isequal(x, t), reversal_equation_vars) # remove iv
    foreach(x -> isparameter(x) && push!(ps, x), reversal_equation_vars)
    return nothing
end

# TODO: turn this into constructor for a `CurrentSystem`
function generate_currents!(gen, dynamics::HodgkinHuxley, synapses, arbor, Vₘ, aₘ)

    currents = Set()
    (; eqs, dvs) = gen

    (; channels, channel_reversals) = dynamics
    paired_channels = zip(channels,
                          map(Base.Fix2(find_reversal,channel_reversals), channels))
    paired_synapses = zip(conductance.(synapses), reversal.(synapses))
    paired_axials = zip(conductances(arbor), reversals(arbor))
    paired_conductances = union(paired_axials, paired_channels, paired_synapses)

    for (chan, Erev) in paired_conductances
        I = IonCurrent(chan)
        g = renamespace(chan, get_output(chan))
        eq = I ~ g * (Vₘ - Erev) * (1 * get_unit(g) isa SpecificConductance ? aₘ : 1)
        validate(eq) && push!(eqs, eq)
        push!(dvs, I)
        push!(currents, I)
    end
    return currents
end


