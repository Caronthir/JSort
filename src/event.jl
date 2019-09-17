struct Scaler
    channel::Int8
    scaler::Int32
    Scaler() = new(0, 0)
    Scaler(x, y) = new(x, y)
end

mutable struct ADCTDC
    channel::Int8
    adc::Int32
    tdc::Int32
    ADCTDC() = new(0, 0, 0)
    ADCTDC(x, y, z) = new(x, y, z)
end

mutable struct ADC
    channel::Int8
    data::Int32
    associated_channel::Int8
    ADC() = new(0, 0, 0)
end


mutable struct Event
    # Channel scaler - Not used
    # scaler::Scaler
    #
    # The SiRi front detector ADC values
    enum::Int8
    e::Array{ADC, 1}

    # The SiRi back detector ADC values
    Δenum::Int8
    Δe::Array{ADC, 1}

    # OSCAR TDC and ADC values
    # Uses a poor man's dict to avoid Julia's costly
    # hashing. labrnum keeps the length of the 
    # labrindices, and labrindices keeps the indices
    # /channels used in labr. Labr has the actual channels
    # and data
    labrnum::Int8 # Keeps track of how many labr are used
    labrindices::Array{Int8, 1} 
    labr::Array{ADCTDC, 1}
    # labr::Dict{UInt32, ADCTDC}

    # Wall clock time
    # time::Union{Int32, Nothing}
    # The TPU pattern
    # pattern::Union{Int32, Nothing}
    function Event()
        e = [ADC() for _ in 1:16]
        Δe = [ADC() for _ in 1:16]
        labr = [ADCTDC() for _ in 1:32]
        labrind = zeros(Int8, 32)
        new(0, e, 0, Δe, 0, labrind, labr)
    end
end

function Base.show(io::IO, e::Event)
    # println("LaBr: $(e.labr)")
    println("# of back: ", e.enum)
    println("# of front: ", e.Δenum)
    println("# of LaBr: ", e.labrnum)
    println("SiRi back(e):")
    [println(adc) for adc in e.e]
    println("\nSiRi front(Δe):")
    [println(adc) for adc in e.Δe]
    println("\nLaBr: ")
    for i in 1:e.labrnum
        chn = e.labrindices[i]
        println("chn: ", chn, " = ", e.labr[chn])
    end
    # println("\nWall clock time: $(repr(e.time))")
    # println("TPU Pattern: $(repr(e.pattern))")
end

function Base.show(io::IO, adc::ADC)
    print("\tchl: $(adc.channel), adc: $(adc.data)")
end

function Base.show(io::IO, adc::ADCTDC)
    print("\tadc: $(adc.adc), tdc: $(adc.tdc)")
end

function setadc!(event::Event, channel, adc)::Nothing
    event.labr[channel+1].adc = adc
    found = false
    for i ∈ 1:event.labrnum
        if event.labrindices[i] == channel+1
            found = true
            break
        end
    end
    if !found
        event.labrnum += 1
        event.labrindices[event.labrnum] = channel+1
        event.labr[channel+1].channel = channel
    end
    nothing
end

function settdc!(event::Event, channel, tdc)::Nothing
    event.labr[channel+1].tdc = tdc
    found = false
    for i ∈ 1:event.labrnum
        if event.labrindices[i] == channel+1
            found = true
            break
        end
    end
    if !found
        event.labrnum += 1
        event.labrindices[event.labrnum] = channel+1
        event.labr[channel+1].channel = channel
    end
    nothing
end

function adde!(event, channel, data)
    event.enum += 1 # First to account for 1-indexing
    event.e[event.enum].channel = channel >> 1
    event.e[event.enum].data = data
end

function addΔe!(event, channel, data)
    event.Δenum += 1
    event.Δe[event.Δenum].channel = channel
    event.Δe[event.Δenum].associated_channel = div(channel, 8)
    event.Δe[event.Δenum].data = data
end

function reset!(e::Event)
    #e.labr = Dict{UInt32, ADCTDC}()
    # Can skip this step, as the 
    # counters decide what is being read
    # for i in 1:16
        # e.e[i].data = e.e[i].channel = 0
        # e.Δe[i].data = e.Δe[i].channel = 0
    # end
    e.enum = 0
    e.Δenum = 0
    for channel in e.labrindices[1:e.labrnum]
        e.labr[channel].adc = 0
        e.labr[channel].tdc = 0
    end
    e.labrnum = 0
end

function eachlabr(e::Event)
    # Returns channel and labr in that channel
    return [(e.labrindices[i], e.labr[e.labrindices[i]]) for i in 1:e.labrnum]
end

function islabrgood(adctdc::ADCTDC)
    return adctdc.adc > 0 &&  adctdc.tdc > 0
end

struct MiniEvent
    front::Int8
    back::Int8
    frontchannel::Int8
    backchannel::Int8
    Δe::Float32
    e::Float32
end

function Base.show(io::IO, event::MiniEvent)
    println("b$(event.back)f$(event.front): Δe = $(event.Δe), e = $(event.e)")
end

struct LaBrEvent
    mevent::MiniEvent
    labr::Array{ADCTDC, 1}
    function LaBrEvent(event::Event, mevent::MiniEvent)
        new(mevent, filter(x -> !islabrgood(x), event.labr[event.labrindices[1:event.labrnum]]))
    end
    function LaBrEvent(mevent::MiniEvent)
        new(mevent, ADCTDC[])
    end
end

function serialize(event::LaBrEvent)
    return [serialize(event.mevent)..., convert(UInt32, length(event.labr))], [serialize(labr) for labr in event.labr]
end

function serialize(event::MiniEvent)
    return (event.front, event.back, event.Δe, event.e) .|> x -> convert(UInt32, x)
end

function serialize(adc::ADCTDC)
    return (adc.channel, adc.adc, adc.tdc) .|> x -> convert(UInt32, x)
end

struct TrueEvent
    f::Int8
    b::Int8
    Δe::Float32
    e::Float32
    labr::Array{ADCTDC, 1}
    function TrueEvent(f, b, de, e, labr)
        new(f, b, de, e, labr)
    end
end
