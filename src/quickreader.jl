using .JSort

using .JSort
import Base: read, read!, fill!, write, close

mutable struct QuickReader <: AbstractReader
    parameters::Parameters
    input::CBuffer
    labr_size::Int
    initialtime::Int
    eγ::Vector{Vector{Float32}}
    source::String
end

function QuickReader(parameters::Parameters)
    j(x) = joinpath(parameters.savepath, x)
    source = joinpath(parameters.readpath, "raw.bin")
    if ~isfile(source)
        error("The file $source does not exist")
    end
    buffer = CBuffer(source)
    egs = [Float32[] for i in 1:33]

    QuickReader(parameters, buffer, 0, -1, egs, source)
end

read(reader::QuickReader, gate::Gate) = read(reader, gate=gate)
read(reader::QuickReader, gate::Gate2D) = read(reader, pgate=gate)

function read(reader::QuickReader; gate::Gate=Gate(), pgate=nothing)
    close(reader.input)
    reader.input = CBuffer(reader.source)

    qkinzpoly = reader.parameters.qkinz
    timecalib = reader.parameters.calibrator[:t]
    timecorrection = reader.parameters.calibrator[:tc]
    gammacalib = reader.parameters.calibrator[:γ]
    ecalib = reader.parameters.calibrator[:e]
    Δecalib = reader.parameters.calibrator[:Δe]
    prompt = reader.parameters.gates[:prompt]
    background = reader.parameters.gates[:background]
    event = Event()
    read_events = front = back = 0
    # TODO Same type as calibrators
    # TODO Buffer2D are not overwritten
    e = Δe = eΔe = Float32(0.0)
    eₓᵗ = eᵧ = Float32(0.0)
    t = tc = Float32(0.0)
    reader.eγ = [Float32[] for i in 1:33]
    for eγ in reader.eγ
        sizehint!(eγ, 2^13)
    end

    while isgood(reader)
        read_events += 1
        fill!(event, reader)
        front, back = event.front, event.back
        e = ecalib[front, back](event.e)
        Δe = Δecalib[front, back](event.de)
        if !isnothing(pgate)
            !pgate(e, Δe) && continue
        end
        # Theoretical excitation energy from kinematics
        eΔe = e + Δe
        # This is a `constant` poly. No calibration needed
        eₓᵗ = qkinzpoly[front](eΔe)
        !gate(eₓᵗ) && continue

        @views for labr in event.labr[1:reader.labr_size]
            labr.channel ∈ reader.parameters.badchannels && continue
            # Gamma detectors
            eᵧ = gammacalib[labr.channel](labr.adc)

            # Time
            if reader.initialtime == -1
                reader.initialtime = event.labr[1].tdc
            end
            t = timecalib[labr.channel](labr.tdc)
            !prompt(t) && continue
            #tc = timecorrection[labr.channel](t, eΔe)
            # Also correct using eᵧ?
            push!(reader.eγ[labr.channel], eᵧ)
        end
    end
end

function close(reader::QuickReader)
    close(reader.input)
end

function ChannelDetector(reader::QuickReader)
    data = GammaDetector[GammaDetector(reader.eγ[chn], chn) for chn in 1:33]
    ChannelDetector(data)
end

plot(r::QuickReader; kwargs...) = plot(ChannelDetector(r); kwargs...)
