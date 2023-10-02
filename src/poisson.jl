struct PoissonDynamics <: AbstractDynamics
    lambda::Num
    function PoissonDynamics(lambda = 0.5)
        @parameters λ = lambda
        return new(λ)
    end
end

poisson_draw(lambda) = rand(Poisson(lambda))

function CompartmentSystem(dynamics::PoissonDynamics; name)
    lambda = dynamics.lambda
    @variables S = 0
    gen = GeneratedCollections(dvs = Set((S,)),
                               ps = Set((lambda,)),
                               eqs = [D(S) ~ 0])

    (; eqs, dvs, ps, observed, systems, defs) = gen
    return CompartmentSystem(nothing, dynamics, nothing, eqs, t, collect(dvs), collect(ps), observed, name,
                             systems, defs, [], Synapse[], Arborization(), StimulusModel[], Point())
end

# TODO: this is kind of a hack
function CompartmentSystem(voltage, dynamics::PoissonDynamics, synapses, arbor, capacitance,
                           geometry, stimuli, defaults, extensions, name)
    return CompartmentSystem(dynamics; name)
end

function Base.convert(::Type{ODESystem}, compartment::CompartmentSystem{PoissonDynamics})
    dvs = get_states(compartment)
    ps = get_ps(compartment)
    eqs = get_eqs(compartment)
    defs = get_defaults(compartment)
    obs = get_observed(compartment)
    syss = convert.(ODESystem, get_systems(compartment))
    return ODESystem(eqs, t, dvs, ps; systems = syss, defaults = defs,
                     name = nameof(compartment), observed = obs, checks = false)
end


