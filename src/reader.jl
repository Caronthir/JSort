using .JSort
import Base: read, read!, fill!, write, close

mutable struct ADCTDC
    channel::Int
    adc::Int
    tdc::Int
    ADCTDC() = new(0, 0, 0)
end

mutable struct Event
    e::UInt32
    de::UInt32
    front::Int
    back::Int
    labr::Vector{ADCTDC}
end

function Event()
    Event(0, 0, 0, 0, ADCTDC[ADCTDC() for i in 1:32])
end

abstract type AbstractReader end

mutable struct Reader <: AbstractReader
    parameters::Parameters
    input::CBuffer
    output::Matrix{CBuffer2D}
    labr_size::Int
    initialtime::Int
    hasread::Bool
    # Syntax is H[istogram] H2[D histogram], followed by variables filling them
    Heγ::  Vector{CBuffer}
    Ht::   Vector{CBuffer}
    Htc::  Vector{CBuffer}
    H2eγt::Vector{CBuffer2D}
    H2eγt_all::   CBuffer2D
    H2eγtc::      CBuffer2D
    H2et::Vector{CBuffer2D}
    H2et_all::   CBuffer2D
    H2etc::      CBuffer2D
    Hexct::Vector{CBuffer}
end

function Reader(parameters::Parameters)
    j(x) = joinpath(parameters.savepath, x)
    output = [CBuffer2D(j("edeb$(b)f$(f).bin"))
               for f in 1:8, b in 1:8]
    fsum   = [CBuffer2D(j("edef$f.bin"))
              for f in 1:8]
    output = [output fsum]
    source = joinpath(parameters.readpath, "raw.bin")
    if ~isfile(source)
        error("The file $source does not exist")
    end
    buffer    = CBuffer(source)

    Heγ       = [CBuffer(j("eg_chn$(chn).bin"), :write) for chn in 1:33]
    Ht        = [CBuffer(j("t_chn$(chn).bin"), :write)  for chn in 1:33]
    Htc       = [CBuffer(j("tc_chn$(chn).bin"), :write) for chn in 1:33]
    H2eγt     = [CBuffer2D(j("eg_t_chn$(chn).bin"))     for chn in 1:33]
    H2eγt_all = CBuffer2D(j("eg_t.bin"))
    H2egtc    = CBuffer2D(j("eg_tc.bin"))
    H2et      = [CBuffer2D(j("e_t_f$(f).bin")) for f in 1:8]
    H2et_all  = CBuffer2D(j("e_t.bin"))
    H2etc     = CBuffer2D(j("e_tc.bin"))
    Hexct     = [CBuffer(j("ex_f$(f).bin"), :write) for f in 1:8]

    Reader(parameters, buffer, output, 0, -1, false,
           Heγ, Ht, Htc, H2eγt, H2eγt_all, H2egtc,
           H2et, H2et_all, H2etc, Hexct)
end

function read(reader::Reader)
    if reader.hasread
        @error "Has already been read"
        return
    end
    qkinzpoly = reader.parameters.qkinz
    timecalib = reader.parameters.calibrator[:t]
    timecorrection = reader.parameters.calibrator[:tc]
    gammacalib = reader.parameters.calibrator[:γ]
    ecalib = reader.parameters.calibrator[:e]
    Δecalib = reader.parameters.calibrator[:Δe]
    egate = reader.parameters.gates[:ex]
    prompt = reader.parameters.gates[:prompt]
    background = reader.parameters.gates[:background]

    event = Event()
    read_events = front = back = 0
    # TODO Same type as calibrators
    # TODO Buffer2D are not overwritten
    e = Δe = eΔe = Float32(0.0)
    eₓᵗ = eᵧ = Float32(0.0)
    t = tc = Float32(0.0)
    while isgood(reader)
        read_events += 1
        fill!(event, reader)
        front, back = event.front, event.back
        # Each individual detector is not calibrated
        write!(reader.output[front, back], event.e, event.de)
        e = ecalib[front, back](event.e)
        Δe = Δecalib[front, back](event.de)

        # Sum over all back for a given front
        write!(reader.output[front, 9], e, Δe)

        # Theoretical excitation energy from kinematics
        eΔe = e + Δe
        # This is a `constant` poly. No calibration needed
        eₓᵗ = qkinzpoly[front](eΔe)

        !egate(eₓᵗ) && continue
        write!(reader.Hexct[front], eₓᵗ )

        @views for labr in event.labr[1:reader.labr_size]
            labr.channel ∈ reader.parameters.badchannels && continue

            # Gamma detectors
            eᵧ = gammacalib[labr.channel](labr.adc)

            # Time
            if reader.initialtime == -1
                reader.initialtime = event.labr[1].tdc
            end
            t = timecalib[labr.channel](labr.tdc)
            tc = timecorrection[labr.channel](t, eΔe)
            # Also correct using eᵧ?

            write!(reader.Heγ[labr.channel], eᵧ)
            write!(reader.Ht[labr.channel], t)
            write!(reader.Htc[labr.channel], tc)
            write!(reader.H2eγt[labr.channel], eᵧ, t)
            write!(reader.H2eγt_all, eᵧ, t)
            write!(reader.H2eγtc, eᵧ, tc)
            write!(reader.H2et[front], e, t)
            write!(reader.H2et_all, e, t)
            write!(reader.H2etc, e, tc)
        end
    end
    @show read_events
    close(reader)
    reader.hasread = true
end

function fill!(event::Event, reader::AbstractReader)
    event.e          = next!(reader)
    event.de         = next!(reader)
    decodeheader!(event, reader, next!(reader))
    getbody!(reader, event.labr)
end

function getbody!(reader::AbstractReader, body::AbstractVector)
    for i in 1:reader.labr_size
        body[i].channel = next!(reader) + 1
        body[i].adc     = next!(reader)
        body[i].tdc     = next!(reader) # Original code has /8
    end
end

function next!(reader::AbstractReader)::UInt32
    next!(reader.input)
end

isgood(reader::AbstractReader) = isgood(reader.input)

function close(reader::Reader)
    map(close, reader.output)
    map(close, reader.Heγ)
    map(close, reader.Ht)
    map(close, reader.H2eγt)
    close(reader.H2eγt_all)
    close(reader.H2eγtc)
    map(close, reader.H2et)
    close(reader.H2et_all)
    close(reader.H2etc)
    close(reader.input)
    map(close, reader.Hexct)
end

function decodeheader!(event::Event, reader::AbstractReader, word::UInt32)
    labrmask  = UInt32(65535)      # Bits 0-15
    backmask  = UInt32(16711680)   # Bits 16-23
    frontmask = UInt32(4278190080) # Bits 23-31
    event.front   = (word & frontmask) >> 24
    event.back    = (word & backmask)  >> 16
    reader.labr_size = word & labrmask
end

