abstract type AbstractDynamics end

struct HodgkinHuxley <: AbstractDynamics
    "Ionic conductances."
    channels::Vector{AbstractConductanceSystem}
    "Equilibrium potentials belonging to ionic membrane conductances."
    channel_reversals::Vector{Num}
end

function filter_reversal!(gen, Erev)
    (; dvs, ps, eqs) = gen
    reversal_equation_vars = Set{Num}()

    if isparameter(Erev)
        push!(ps, Erev)
    else
        push!(dvs, Erev)
        # FIXME: hacked solution (assumes eq for dynamic reversal stored in default value)
        if MTK.hasdefault(Erev)
            rhs = getdefault(Erev)
            get_variables!(reversal_equation_vars, rhs)
            push!(eqs, Erev ~ rhs)
        end
    end
    setdiff!(reversal_equation_vars, t) # remove iv
    foreach(x -> isparameter(x) && push!(ps, x), reversal_equation_vars)
    return
end

# FIXME!!!!
function generate_currents!(currents, gen, paired_conductances, Vₘ, aₘ)
    CurrentSystem()
#    (; eqs, dvs) = gen
#    for (chan, Erev) in paired_conductances
#        I = IonCurrent(chan)
#        g = renamespace(chan, get_output(chan))
#        eq = I ~ g * (Vₘ - Erev) * (1 * get_unit(g) isa SpecificConductance ? aₘ : 1)
#        validate(eq) && push!(eqs, eq)
#        push!(dvs, I)
#        push!(currents, I)
#    end
    return
end

function generate_currents!(currents, gen, dynamics::HodgkinHuxley, Vₘ, aₘ)
    (; channels, channel_reversals) = dynamics
    paired_channels = zip(channels,
                          map(Base.Fix2(find_reversal,channel_reversals), channels))
    generate_currents!(currents, gen, paired_conductances, Vₘ, aₘ)
end

function generate_currents!(currents, gen, synapses::Vector{Synapse}, Vₘ, aₘ)
    paired_synapses = zip(conductance.(synapses), reversal.(synapses))
    generate_currents!(currents, gen, paired_synapses, Vₘ, aₘ)
end
       
function generate_currents!(currents, gen, arbor::Arborization, Vₘ, aₘ)
    paired_axials = zip(conductances(arbor), reversals(arbor))
    generate_currents!(currents, gen, paired_axials, Vₘ, aₘ)
end
