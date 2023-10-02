struct PoissonDynamics <: AbstractDynamics
    lambda::Num
    function PoissonDynamics(lambda = 0.5)
        @parameters λ = lambda
        return new(λ)
    end
end

poisson_draw(lambda) = rand(Poisson(lambda))
Symbolics.@register_symbolic poisson_draw(x)

function CompartmentSystem(dynamics::PoissonDynamics; name)
    lambda = dynamics.lambda
    @variables S = 0
    gen = GeneratedCollections(dvs = Set((S,)),
                               ps = Set((lambda,)),
                               eqs = [S ~ poisson_draw(lambda)])

    (; eqs, dvs, ps, observed, systems, defs) = gen
    return CompartmentSystem(S, dynamics, nothing, eqs, t, collect(dvs), collect(ps), observed, name,
                             systems, defs, [], Synapse[], Arborization(), StimulusModel[], Point())
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


