
abstract type Stimulus end
get_stimulus(x::Num) = getmetadata(x, Stimulus)

# Fixed value
struct Bias{T} <: Stimulus
    val::T
end

function Bias(amplitude::T; name = Base.gensym("bias")) where {T <: Current}
    val = ustrip(µA, amplitude)
    sym = only(@parameters $name=val [unit = µA])
    return setmetadata(sym, Stimulus, Bias{T}(amplitude))
end

struct PulseTrain{U} <: Stimulus
    interval::Any
    duration::Any
    delay::Any
    npulses::Int
    offset::Any
    amplitude::Any
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
    val = ustrip(µA, offset)
    sym = only(@variables $name(t)=val [unit = µA])
    return setmetadata(sym, Stimulus,
                       PulseTrain{T1}(ustrip(ms, interval),
                                      ustrip(ms, duration),
                                      ustrip(ms, delay),
                                      npulses,
                                      val, # offset
                                      ustrip(µA, amplitude)))
end

# Arbitrary input
#struct CommandCurrent end
#struct CommandPotential end
#struct Stimulator{T<:Stimulus}
#    channels::Vector{T} 
#end
#
#struct Waveform{T} <: Stimulus end
