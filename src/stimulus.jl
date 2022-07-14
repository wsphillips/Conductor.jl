
abstract type AbstractStimulus end

struct Stimulator{T<:AbstractStimulus}
    channels::Vector{T} 
end

struct PulseTrain{T<:Real} <: AbstractStimulus
    interval::T
    duration::T
    delay::T
    npulses::Int
    offset::T
    amplitude::T
end

struct Waveform{T<:Real} <: AbstractStimulus end

