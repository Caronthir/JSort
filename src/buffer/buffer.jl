const lib = "./libbuffer"
abstract type Buffer end
mutable struct Buffer2D <: Buffer end

function makebuffer(path::AbstractString, mode::Symbol=:read)
    ptr = ccall((:makeBufferR, lib), Ptr{Buffer}, (Cstring, ), path)
end

function close(buffer::Ptr{Buffer})
    ccall((:close, lib), Cvoid, (Ptr{Buffer},), buffer)
end

function next!(buffer::Ptr{Buffer})
    val = ccall((:next, lib), UInt32, (Ptr{Buffer},), buffer)
end

function isgood(buffer::Ptr{Buffer})
    good = ccall((:good, lib), Bool, (Ptr{Buffer},), buffer)
end

function makebuffer2d(path::AbstractString, mode::Symbol=:write)
    ptr = ccall((:makeBufferW2D, lib), Ptr{Buffer2D}, (Cstring, ), path)
end

function write!(buffer::Ptr{Buffer2D}, mat::AbstractMatrix)
    for row in 1:size(mat, 1)
        #println("Writing [$(mat[row, 1]), $(mat[row, 2])]")
        ccall((:fill2D, lib), Cvoid, (Ptr{Buffer2D}, UInt32, UInt32),
              buffer, mat[row, 1], mat[row, 2])
    end
end

function close(buffer::Ptr{Buffer2D})
    ccall((:close2D, lib), Cvoid, (Ptr{Buffer2D},), buffer)
end

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
    @time mat = [r+c for r in 1:100_000_000, c in 0:1]
    @time write!(buffer, mat)
    @time close(buffer)
end

@time main3()
