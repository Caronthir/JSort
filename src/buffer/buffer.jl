#module CBuffers
#export AbstractCBuffer, CBuffer, CBuffer2D, open, close, isgood, write!, next!
import Base:open, close
using Base.Filesystem: touch, abspath

#TODO FIX!
#const lib = abspath("./libbuffer")
const lib = "$(homedir())/gits/JSort/src/buffer/libbuffer.so"
abstract type Buffer end
struct Buffer2D <: Buffer end

function makebuffer(path::AbstractString; mode::Symbol=:read)
    if mode == :read
        ptr = ccall((:makeBufferR, lib), Ptr{Buffer}, (Cstring, ), path)
    elseif mode == :write
        ptr = ccall((:makeBufferW, lib), Ptr{Buffer}, (Cstring, ), path)
    else
        error("Only :read and :write supported for CBuffer")
    end
end

function close(buffer::Ptr{Buffer})
    ccall((:close, lib), Cvoid, (Ptr{Buffer},), buffer)
end

function next!(buffer::Ptr{Buffer})::UInt32
    val = ccall((:next, lib), UInt32, (Ptr{Buffer},), buffer)
end

function isgood(buffer::Ptr{Buffer})
    good = ccall((:good, lib), Bool, (Ptr{Buffer},), buffer)
end

function isgood(buffer::Ptr{Buffer2D})
    good = ccall((:good, lib), Bool, (Ptr{Buffer2D},), buffer)
end

function makebuffer2d(path::AbstractString, mode::Symbol=:write)
    ptr = ccall((:makeBufferW2D, lib), Ptr{Buffer2D}, (Cstring, ), path)
end

function write!(buffer::Ptr{Buffer}, x::UInt32)
    ccall((:fill, lib), Cvoid, (Ptr{Buffer}, UInt32), buffer, x)
end
write!(buffer::Ptr{Buffer}, x) = write!(buffer, reinterpret(UInt32, x))

function write!(buffer::Ptr{Buffer2D}, mat::AbstractMatrix)
    for row in 1:size(mat, 1)
        #println("Writing [$(mat[row, 1]), $(mat[row, 2])]")
        ccall((:fill2D, lib), Cvoid, (Ptr{Buffer2D}, UInt32, UInt32),
              buffer, mat[row, 1], mat[row, 2])
    end
end

function write!(buffer::Ptr{Buffer2D}, x::UInt32, y::UInt32)
    ccall((:fill2D, lib), Cvoid, (Ptr{Buffer2D}, UInt32, UInt32),
          buffer, x, y)
end
write!(buffer::Ptr{Buffer2D}, x, y) = write!(buffer, reinterpret(UInt32, x),
                                             reinterpret(UInt32, y))

function close(buffer::Ptr{Buffer2D})
    ccall((:close2D, lib), Cvoid, (Ptr{Buffer2D},), buffer)
end

abstract type AbstractCBuffer end
struct CBuffer <: AbstractCBuffer
    ptr::Ptr{Buffer}
    mode::Symbol
end
struct CBuffer2D <: AbstractCBuffer
    ptr::Ptr{Buffer2D}
    mode::Symbol
end
CBuffer(path::AbstractString, mode::Symbol=:read) = open(path, CBuffer, mode=mode)
CBuffer2D(path::AbstractString) = open(path, CBuffer2D)

function open(path::AbstractString, T::Type{AbstractCBuffer}; mode::Symbol=:read)
    open(path, T, mode=Val{mode})
end

function open(path::AbstractString, T::Type{CBuffer}; mode=:read)
    if mode != :read && mode != :write
        error("Only reading xor writing supported for CBuffer")
    end
    buffer_ptr = makebuffer(path, mode=mode)
    return CBuffer(buffer_ptr, mode)
end

function open(path::AbstractString, T::Type{CBuffer2D}; mode::Symbol=:write)
    if mode != :write
        error("Only writing supported for CBuffer2D")
    end
    touch(path)  # The file is not created by the C-code
    buffer_ptr = makebuffer2d(path)
    return CBuffer2D(buffer_ptr, :write)
end

function close(buffer::AbstractCBuffer)
    close(buffer.ptr)
end

function next!(buffer::AbstractCBuffer)
    next!(buffer.ptr)
end

write!(buffer::CBuffer2D, m::AbstractMatrix) = write!(buffer.ptr, m)
write!(buffer::CBuffer2D, x, y) = write!(buffer.ptr, x, y)
write!(buffer::CBuffer, x) = write!(buffer.ptr, x)
isgood(buffer::AbstractCBuffer) = isgood(buffer.ptr)
isgood(buffer::CBuffer) = isgood(buffer.ptr)

function main()
    buffer = makebuffer("duck")
    fact = read("duck") |> x -> reinterpret(UInt32, x)
    i = 1
    test = similar(fact)
    @allocated while isgood(buffer)
        test[i] = next!(buffer)
        if test[i] != fact[i]
            println(i)
        end
        i += 1
    end
    close(buffer)
end

function main2()
    buffer = makebuffer2d("ede.bin")
    mat = UInt32[
        10 2
        11 3
        4 5
        12 6]
    write!(buffer, mat)
    close(buffer)
end

function main3()
    buffer = makebuffer2d("ede.bin")
    @time mat = [r+c for r in 1:10_00_000, c in 0:1]
    @time write!(buffer, mat)
    @time close(buffer)
end

#end
