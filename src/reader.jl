using .JSort
using ArgCheck
using Base.Iterators
using Distributed
using ProgressMeter
using Printf
using Profile
using Traceur
import Base.read

function processfile(path::T, parameters::Parameters) where T<:AbstractString
    if nworkers() > 1
        processfile_parallel(path, parameters)
    else
        processfile_serial(path, parameters)
    end
end

function processfile_parallel(path, parameters::Parameters)
    filesize = stat(path).size
    contentsize = ceil(Int64, parameters.percent/100*filesize)

    # Distribute the chunks of the file to read
    partitions = partitionfile(path, contentsize)
    progresschannel = RemoteChannel(()->Channel{Int64}(256))
    @async showprogress(progresschannel, length(workers())*10)
    print("Constructing matrices...")
    matrices = creatematrices(SharedArray)
    println("done.")

    responses = Vector{Any}(undef, nworkers())
    @sync for (i, worker) in enumerate(workers())
        @async responses[i] = remotecall_fetch(processfilechunk!, worker, matrices,
                               path, partitions[i], partitions[i+1],
                               progresschannel, parameters)
    end
    N = 0
    S, U, R = 0.0, 0.0, 0.0
    for (n, s, u, r) in responses
        # @show (n, s, u, r)
        N += n
        S += s/length(responses) * 100
        U += u/length(responses) * 100
        R += r/length(responses) * 100
    end
    @printf "Managed to read %g events. %.1f %% of the data was usable, where %.1f %% was unreadable and %.1f %% rejected before sorting" N S U R
    return matrices
end

function processfile_serial(path::AbstractString, parameters::Parameters)
    filesize = stat(path).size
    contentsize = ceil(Int64, parameters.percent/100*filesize)
    progresschannel = Channel{Int64}(256)
    # AARGH! Why doesn't progressmeter work!?
    # Changed the size to 2 to force the thread to write
    @async showprogress(progresschannel, 10)
    labrevents, numevents, successrate, unreadable, rejected = processfilechunk!(path, 0, contentsize, progresschannel, parameters)
    @printf "Managed to read %g events. %.1f %% of the data was usable, where %.1f %% was unreadable and %.1f %% rejected before sorting" numevents 100successrate 100unreadable 100rejected
    return labrevents
end

function processfilechunk!(path, start, stop, progresschannel, parameters::Parameters)
    #=
    Algorithm: Go the the chunk in the file and read the bytes as unsigned int 32.
    Step through the words to find header words. The header contains the length
    of the data packet. Send the data packet to the unpacker which decodes
    the entire packet. If the data unpacking fails, or if the event is rejected by some
    simple heuristics, the rest of the packet is ignored. Otherwise the decoded packet
    is sent as an event to the sorter. 
    Since the data can be corrupt, the corrupt portions will be ignored until a new
    header is found.
    =#

    fhandle = open(path)
    seek(fhandle, start)
    # Make sure the conversion from UInt8 to UInt32 is possible
    readsize = stop-start - mod(stop-start, 4)
    words = read(fhandle, readsize) |> x -> reinterpret(UInt32, x)
    # println("Reading $start to $stop and $(length(words)) words")

    # Counters for summarizing the quality of the data
    numevents = 0
    validwords = 0
    unreadable_words = 0
    rejected_words = 0

    # Progress marker
    prev = length(words)/10

    i = 1
    event = Event()
    event_counter = 1
    labrevents = LaBrEvent[]
    while i < length(words)
        if isheader(words[i])
            packetsize = ndw(words[i])
            if i+packetsize > length(words)
                # println("\nIncomplete chunk. Aborting")
                break
            end

            success = unpack!(event, words[i+1:i+packetsize])
            if success
                #sortevent!(matrices, event, parameters)
                labr = correlateevent(event)
                if !isnothing(labr)
                    push!(labrevents, labr)
                    numevents += 1
                    validwords += packetsize+1
                end
            else
                rejected_words += packetsize
            end
            i += packetsize+1
        elseif isendofbuffer(words[i])
            i+=1
            validwords += 1
        else
            # Something is wrong. Ignore rest of packet
            # This seems to always occur after an EOB and is of constant length
            while i < length(words) && !isheader(words[i])
                i+=1
                unreadable_words += 1
            end
        end

        if i > prev
            # Update at each 10% milepæl
            put!(progresschannel, i)
            prev += floor(length(words)/10)
        end
    end
    # Mark the progress as done
    while i < length(words)
        put!(progresschannel, i)
        i += floor(length(words)/10)
    end
    close(fhandle)
    return labrevents, numevents, validwords/length(words), unreadable_words/length(words), rejected_words/length(words)
end

function unpack!(event::Event, packet)::Bool
    reset!(event)
    for i in eachindex(packet)
        word = packet[i]

        # Be wary of corrupt data
        if boe(word) ≠ 0x00
            # return nothing
            continue
        end

        # By stupid data design, the next word is required to
        # determine the quality of the current word
        nextword::Union{Missing, UInt32} = if i+1 < length(packet) packet[i+1] else missing end
        if(!decodeword!(event, word, nextword))
            return false 
        end
    end

    # Reject the event if there are 2 ΔE events in the same detector
    # i.e pileup
    if event.enum > 0 && event.Δenum > 1
        removepileup!(event)
    end

    # Reject the event if either the front teller or back teller is empty
    if event.enum == 0 || event.Δenum == 0
        return false
    end
    return true
end

# Extract the header
boe(x::UInt32) = (x&(0xC0000000::UInt32)) >> 28
# The length of the packet given in the header
ndw(x::UInt32) =  x&(0x000000ff::UInt32)
# Which box has fired
box(x::UInt32) = (x&(0x3f800000::UInt32)) >> 23
# The channel
chn(x::UInt32) = (x&(0x007f0000::UInt32)) >> 16
# The actual data
dta(x::UInt32) =  x&(0x0000ffff::UInt32)
# Check if is header
isheader(x::UInt32) = boe(x) == 0xC
# Check if is end of buffer
isendofbuffer(x::UInt32) = x == 0x80000000
# Check if is a guard ring on the back counter
isguardring(channel::UInt32) = (channel&1) == 0 || channel ≥ 16

function decodeword!(event::Event, word, nextword)::Bool
    boxnum = box(word)
    channel = chn(word)
    data = dta(word)

    if boxnum == 0x00
        # TPU pattern ch 0 - 3
        # if channel == 0 && isnothing(event.pattern)
            # event.pattern = data
        # else
            # return true
        # end
    elseif boxnum == 0x01
        # Wall-clock time ch 16 (high) and 17 (low)
        # This continuously fails
        # if channel ≠ 16 || ismissing(nextword) || chn(nextword) ≠ 17
            # return true
        # end
        # time = data << 16 | dta(nextword::UInt32)
        # event.time = time
    elseif boxnum == 0x02
        # if ismissing(nextword)
            # return true
        # end
        # # VME scaler 1151N ch 0-15 (low) and 16-31 (hi)
        # # Doesn't seem to be used
        # return true
        # if chn(nextword) - channel ≠ 16
            # return true
        # end
        # scaler = Scaler(channel, dta(nextword) << 16 | data)
        # if event.scaler.scaler != 0
            # error("Fix array")
        # end
        # event.scaler = scaler
    elseif boxnum == 0x10
        # Time of NaI ch 0-31
        settdc!(event, channel, data)
    elseif boxnum == 0x20 || boxnum == 0x24
        # Energy of LaBr ch 0-31, MADC ch 0-31
        # Reject if adc < 0
        data < 0 && return true
        setadc!(event, channel, data)
    elseif boxnum == 0x21
        # Energy of E ch 0-31
        isguardring(channel) && return true
        adde!(event, channel, data)
    elseif boxnum == 0x22
        # Energy of ΔE1 ch 0-31
        # try
          addΔe!(event, channel, data)
    #     catch e
    #         @show e
    #         @show event
    #     end
    # elseif boxnum == 0x23
        # Energy of ΔE2 ch 32-61
        addΔe!(event, channel+32, data)
    else
        return false
    end
    return true
end

function removepileup!(event::Event)
    # If an event has multiple Δe corresponding to the same back detector, remove them

    # Find all duplicates and store their index
    bad_id = Set{Int8}()
    final = event.Δenum
    for i in Int8(1):Int8(final)
        for j in Int8(i+1):Int8(final)
            if event.Δe[i].associated_channel == event.Δe[j].associated_channel
                push!(bad_id, i, j)
            end
        end
    end
    if length(bad_id) == 0
        return
    elseif length(bad_id) == final
        event.Δenum = 0
        return
    end
    # Move all elements down, except those with
    # bad index. A byproduct is that the
    # elements not overwritten remains,
    # but this should not interfere as the
    # counter shows the correct end
    event.Δenum = 0
    for i in Int8(1):Int8(final)
        if i ∉ bad_id
            event.Δenum += 1
            if i ≠ event.Δenum
                event.Δe[event.Δenum].channel = event.Δe[i].channel
                event.Δe[event.Δenum].data = event.Δe[i].data
                event.Δe[event.Δenum].associated_channel = event.Δe[i].associated_channel
            end
        end
    end
end


function correlate(event::Event)::Union{Nothing, MiniEvent}
    foundevent = false
    let minievent
        for Δei in 1:event.Δenum
            Δe = event.Δe[Δei]
            for ei in 1:event.enum
                e = event.e[ei]
                back = e.channel
                if Δe.associated_channel ≠ back
                    continue
                end
                front = Δe.channel % 8
                # We can only handle 1 correlated event
                if foundevent
                    return nothing
                end
                foundevent = true
                minievent = MiniEvent(front+1, back+1, Δe.channel, e.channel, Δe.data, e.data)
            end
        end
        if foundevent
            return minievent
        else
            return nothing
        end
    end
end

function correlateevent(event::Event)
    minievent = correlate(event)
    isnothing(minievent) && return
    labrs = LaBrEvent(minievent)
    isgood = false
    for (channel, labr) in eachlabr(event)
        if islabrgood(labr)
            isgood = true
            push!(labrs.labr, ADCTDC(labr.channel, labr.adc, labr.tdc))
        end
    end
    !isgood && return
    return labrs
end

function showprogress(channel, max)
    progress = Progress(max, "Unpacking events: ")
    count = 0
    while count < max
        try
            take!(channel)
            next!(progress)
        catch
            println("opsie")
            break
        end
        count += 1
    end
end

function partitionfile(path, readsize)
    step = floor(Int64, readsize/nworkers())
    partitions = collect(0:step:readsize+1)

    # Move the endpoints to find the first header
    fhandle = open(path)
    partitions₂ = [partitions[1]]
    for pos in partitions[2:end-1]
        seek(fhandle, pos)
        while !isheader(read(fhandle, UInt32))
            pos -= 1
            seek(fhandle, pos)
        end
        push!(partitions₂, pos)
    end
    push!(partitions₂, partitions[end])
    partitions₂
end

function measurequality(path)
    fhandle = open(path)
    stop = stat(path).size
    seek(fhandle, 0)
    # Make sure the conversion from UInt8 to UInt32 is possible
    readsize = stop - mod(stop, 4)
    words = read(fhandle, readsize) |> x -> reinterpret(UInt32, x)
    # println("Reading $start to $stop and $(length(words)) words")
    wordtype = zeros(Int8, length(words)) .- Int8(1)
    println(typeof(wordtype))

    # Counters for summarizing the quality of the data
    numevents = 0
    validwords = 0
    unreadable_words = 0
    rejected_words = 0

    HEADER = 0
    BODY = 1
    REJECTED = 2
    EOB = 3
    UNREADABLE = 4

    # Progress marker
    prev = length(words)/10
    progress = Progress(10)

    i = 1
    while i < length(words)
        if isheader(words[i])
            wordtype[i] = HEADER
            packetsize = ndw(words[i])
            if i+packetsize > length(words)
                # println("\nIncomplete chunk. Aborting")
                break
            end

            event = unpack(words[i+1:i+packetsize])
            if !isnothing(event)
                numevents += 1
                validwords += packetsize+1
                wordtype[i+1:i+packetsize] .= BODY
            else
                rejected_words += packetsize
                wordtype[i+1:i+packetsize] .= REJECTED
            end
            i += packetsize+1
        elseif isendofbuffer(words[i])
            i+=1
            validwords += 1
            wordtype[i] = EOB
        else
            # Something is wrong. Ignore rest of packet
            while i < length(words) && !isheader(words[i])
                i+=1
                unreadable_words += 1
                wordtype[i] = UNREADABLE
            end
        end

        if i > prev
            # Update at each 10% milepæl
            next!(progress)
            prev += floor(length(words)/10)
        end
    end
    # Mark the progress as done
    while i < length(words)
        next!(progress)
        i += floor(length(words)/10)
    end
    write("quality.dat", wordtype)
end
