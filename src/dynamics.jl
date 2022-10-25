abstract type AbstractDynamics end

struct HodgkinHuxley <: AbstractDynamics
    "Ionic conductances."
    channels::Vector{AbstractConductanceSystem}
    "Equilibrium potentials belonging to ionic membrane conductances."
    channel_reversals::Vector{Num}
end

function reversals(dynamics::HodgkinHuxley)
    (; channels, channel_reversals) = dynamics
    map(Base.Fix2(find_reversal,channel_reversals), channels)
end

conductances(dynamics::HodgkinHuxley) = dynamics.channels

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

function generate_currents!(current_systems, gen, dynamics, Vₘ, aₘ)
    paired_conductances = zip(conductances(dynamics), reversals(dynamics))
    for (cond, Erev) in paired_conductances
        filter_reversal!(gen, Erev)
        sys = CurrentSystem(Vₘ, cond, Erev; aₘ = aₘ)
        push!(current_systems, sys) 
        push!(gen.systems, sys)
    end
    return
end

function generate_currents!(current_systems, gen, stimuli::Vector{<:Stimulus}, Vₘ, aₘ)
    for stimulus in stimuli
        sys = CurrentSystem(Vₘ, stimulus)
        push!(current_systems, sys)
        push!(gen.systems, sys)
    end
end
