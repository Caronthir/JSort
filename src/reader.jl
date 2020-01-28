using .JSort
import Base: read, read!, fill!, write, close


mutable struct ADCTDC
    channel::Int
    adc::Int
    tdc::Int
    ADCTDC() = new(0, 0, 0)
end

mutable struct Event
    e::Int
    de::Int
    front::Int
    back::Int
    labr::Vector{ADCTDC}
end

function Event()
    Event(0, 0, 0, 0, ADCTDC[ADCTDC() for i in 1:32])
end

abstract type Buffer end

mutable struct Buffer2D <: Buffer
    source::IO
    buffer::Matrix{UInt32}
    bp::Int
    read::Int
end
function Buffer2D(io::IO; size=2^16)
    Buffer2D(io, Matrix{UInt32}(undef, size, 2), 1, 0)
end
Buffer2D(path::AbstractString; size=2^16, rw="w") = Buffer2D(open(path, rw), size=size)

mutable struct Buffer1D <: Buffer
    source::IO
    buffer::Vector{UInt32}
    bp::Int
    read::Int
end
function Buffer1D(io::IO; size=2^16)
    Buffer1D(io, Vector{UInt32}(undef, size), 1, 0)
end
Buffer1D(path::AbstractString; size=2^16, rw="w") = Buffer1D(open(path, rw), size=size)

mutable struct Reader
    buffer::Buffer1D
    output::Matrix{Buffer2D}
    labr_size::Int
end

function Reader(parameters::Parameters; buffersize=2^15)
    output = [Buffer2D(joinpath(parameters.savepath, "edeb$(b)f$(f).bin"))
               for f in 1:8, b in 1:8]
    fsum   = [Buffer2D(joinpath(parameters.savepath, "edef$f.bin"))
              for f in 1:8]
    output = [output fsum]
    source = open(joinpath(parameters.readpath, "output.bin"), "r")
    buffer = Buffer1D(source, size=buffersize)
    Reader(buffer, output, 0)
end

function read(reader::Reader)
    event = Event()
    fillbuffer!(reader)
    read_events = 0
    while isgood(reader)
        fill!(event, reader)
        fill!(reader.output[event.front, event.back], event.e, event.de)
        # Sum over all back for a given front
        fill!(reader.output[event.front, 9], event.e, event.de)
        read_events += 1
    end
    @show read_events
    write(reader)
    close(reader)
end

function fill!(event::Event, reader::Reader)
    event.e = next!(reader)
    event.de = next!(reader)
    event.front = next!(reader)
    event.back = next!(reader)
    reader.labr_size = next!(reader)
    getbody!(reader, event.labr)
end

function fillbuffer!(reader::Reader)
    fillbuffer!(reader.buffer)
end

function getheader(reader::Reader)
    r = next!(reader.buffer, 5)
    reader.labr_size = r[5]
    r
end

function getbody!(reader::Reader, body)
    for i in 1:reader.labr_size
        body[i].channel = next!(reader)
        body[i].adc     = next!(reader)
        body[i].tdc     = next!(reader)
    end
end


function fillbuffer!(buffer::Buffer1D)
    buffer.bp = 1
    remaining = stat(buffer.source).size - buffer.read
    # Resize the buffer to the remaining space if there isn't enough
    # data left to fill it
    if remaining >= sizeof(buffer.buffer)
        read!(buffer.source, buffer.buffer)
    else
        println("Last piece")
        remaining = remaining / sizeof(eltype(buffer.buffer))|> Int
        println(remaining)
        buffer.buffer = buffer.buffer[1:remaining]
        read!(buffer.source, buffer.buffer)
    end
    buffer.read += sizeof(buffer.buffer)
end

function isgood(buffer::Buffer)::Bool
    return !eof(buffer.source) || buffer.bp < length(buffer.buffer)-1
end
isgood(r::Reader) = isgood(r.buffer)

function conditionalfill!(buffer::Buffer)
    if buffer.bp <= length(buffer.buffer)
        return
    end
    fillbuffer!(buffer)
end

function next!(buffer::Buffer)
    val = buffer.buffer[buffer.bp]
    buffer.bp += 1
    conditionalfill!(buffer)
    return val
end

function next!(buffer::Buffer, size::Int)
    [next!(buffer) for _ in 1:size]
end
next!(r::Reader) = next!(r.buffer)
next!(r::Reader, i::Int) = next!(r.buffer, i)

function fill!(buffer::Buffer, x, y)
    #println("Filling! $x, $y")
    buffer.buffer[buffer.bp, 1] = x
    buffer.buffer[buffer.bp, 2] = y
    buffer.bp += 1
    conditionalwrite(buffer)
end

function conditionalwrite(buffer::Buffer)
    if buffer.bp == size(buffer.buffer, 1)
        write(buffer)
    elseif buffer.bp > size(buffer.buffer, 1)
        throw("Corrupt buffer: index too large")
    end
end

function write(buffer::Buffer2D)
    if buffer.bp != 1
      @views write(buffer.source, transpose(buffer.buffer[1:buffer.bp, :]))
      buffer.bp = 1
    end
end

function close(buffer::Buffer)
    close(buffer.source)
end

function write(reader::Reader)
    map(write, reader.output)
end

function close(reader::Reader)
    map(close, reader.output)
    close(reader.buffer)
end
