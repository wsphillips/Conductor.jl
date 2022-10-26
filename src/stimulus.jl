
abstract type Stimulus end

struct Bias{T} <: Stimulus
    val::T
    name::Symbol
end

function Bias(amplitude::T; name = Base.gensym("bias")) where {T <: Current}
    return Bias{T}(amplitude, name)
end

struct PulseTrain{U} <: Stimulus
    interval::Any
    duration::Any
    delay::Any
    npulses::Int
    offset::Any
    amplitude::Any
    name::Symbol
end

function current_pulses(t, x)
    (; interval, duration, delay, npulses, offset, amplitude) = x
    if t > delay
        telapsed = t - delay
        pulse_number = (telapsed ÷ interval) + 1
        pulse_number > npulses && return offset
        t_in_period = telapsed - (pulse_number - 1) * interval
        t_in_period < duration && return amplitude
    end
    return offset
end

@register_symbolic current_pulses(t, x::PulseTrain)

function PulseTrain(; amplitude::T1,
                    duration::Time,
                    offset::T2 = 0.0µA,
                    delay::Time,
                    npulses = 1,
                    interval::Time = duration,
                    name = Base.gensym("pulsetrain")) where {T1 <: Current, T2 <: Current}
    return PulseTrain{T1}(ustrip(ms, interval),
                          ustrip(ms, duration),
                          ustrip(ms, delay),
                          npulses,
                          ustrip(µA, offset),
                          ustrip(µA, amplitude),
                          name)
end

function IonCurrent(stim::Bias{T}) where {T <: Current}
    IonCurrent(NonIonic, stim.val*µA; dynamic = false, name = :Iₑ)
end

function IonCurrent(stim::PulseTrain{T}) where {T <: Current}
    IonCurrent(NonIonic, stim.offset*µA; dynamic = true, name = :Iₑ)
end

# Arbitrary input
#struct CommandCurrent end
#struct CommandPotential end
#struct Stimulator{T<:Stimulus}
#    channels::Vector{T} 
#end
#
#struct Waveform{T} <: Stimulus end
