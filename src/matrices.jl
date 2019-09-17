using SharedArrays
using .JSort

const MAX_E = 20_000
const MAX_ΔE = 8_000
const NUMBINS_E = 1000
const NUMBINS_CHANNEL = 2000

struct OMatrix{T}
    matrix::T
    minx::Int64
    miny::Int64
    maxx::Int64
    maxy::Int64
    α::Array{Float64, 1}
    title::String
    xlabel::String
    ylabel::String

    function OMatrix{SharedArray{Int64, 2}}(shape; max=(), min=(0, 0), title="", xlabel="", ylabel="")
        max = if length(max) == 2 max else shape end
        α = collect(shape./(max.-min))
        new(SharedArray{Int64, 2}(shape, init=0), min[1], min[2], max[1], max[2], α, title, xlabel, ylabel)
    end

    function OMatrix{SharedArray{Int64, 1}}(shape::T where T<:Number; min=0, max=-1, title="", xlabel="", ylabel="counts")
        max = if max ≠ -1 max else shape end
        α = [shape/(max-min), 0]
        new(SharedArray{Int64, 1}(shape, init=0),
            min, -1, max, -1, α, 
            title, xlabel, ylabel)
    end

    function OMatrix{Array{Int64, 2}}(shape; max=(), min=(0, 0), title="", xlabel="", ylabel="")
        max = if length(max) == 2 max else shape end
        α = collect(shape./(max.-min))
        new(zeros(shape), min[1], min[2], max[1], max[2], α, title, xlabel, ylabel)
    end

    function OMatrix{Array{Int64, 1}}(shape::T where T<:Number; min=0, max=-1, title="", xlabel="", ylabel="counts")
        max = if max ≠ -1 max else shape end
        α = [shape/(max-min), 0]
        new(zeros(shape),
            min, -1, max, -1, α, 
            title, xlabel, ylabel)
    end
end

function creatematrices(T)#::Dict{Symbol, Union{OMatrix{T{Int64, 1}}, OMatrix{T{Int64, 2}}}}
    # TODO Since OMatrix is parametric, matrices[] is no longer type stable,
    # and functions called on an element of matrices[] must use dynamic dispatch.
    H1 = T{Int64, 1}; H2 = T{Int64, 2}
    matrices = Dict{Symbol, OMatrix}()
    matrices[:back]     = OMatrix{H2}((NUMBINS_CHANNEL, 8), max=(MAX_E, 8), 
                              title="Back detector energies", xlabel="SiRi Energy [keV]",
                              ylabel="Hetector Channel")
    matrices[:front]    = OMatrix{H2}((NUMBINS_CHANNEL, 64), max=(MAX_ΔE, 64),
                               title="Front detector energies", xlabel="SiRi Energy [keV]",
                               ylabel="Hetector Channel")
    matrices[:eΔe]      = OMatrix{H2}((NUMBINS_E, NUMBINS_E), max=(MAX_E, MAX_ΔE),
                               title=raw"Total $E$-$\Helta E$", xlabel=raw"$E$ [keV]",
                               ylabel=raw"$\Helta E$ [keV]")
    matrices[:emul]     = OMatrix{H1}(10)
    matrices[:Δemul]    = OMatrix{H1}(64)
    matrices[:heΔe]     = OMatrix{H1}(1000, max=MAX_E)
    matrices[:thick]    = OMatrix{H1}(500)
    matrices[:eΔethick] = OMatrix{H2}((1000, 1000), max=(MAX_E, MAX_ΔE))
    matrices[:hex]      = OMatrix{H1}(500, max=15_000, min=-500)
    matrices[:labrmul]  = OMatrix{H1}(32, title="LaBr Multplicity", xlabel="Multiplicity")


    # for (b, f) in Iterators.product(1:8, 1:8)
        # matrices[Symbol("eΔeb", b, "f", f)]  = OMatrix((NUMBINS_E, NUMBINS_E), max=(MAX_E, MAX_ΔE))
        # matrices[Symbol("heΔeb", b, "f", f)] = OMatrix(2000, max=MAX_E)
    # end
    # for f in 1:8
        # matrices[Symbol("heΔerf$f")] = OMatrix(1000, max=MAX_E)
        # matrices[Symbol("hexrf$f")]  = OMatrix(500, min=-500, max=15_000)  # NOTE Why 15_000?
    # end
    return matrices
end

function save(dict, path)
    # TODO: Strings are unreadable in Python
    # TODO: Clean up and specialize this function
    !isdir(path) && mkdir(path)
    for (key, matrix) in dict
        name = replace(String(key), "Δ"=>"d")*".jld"
        if name ∉ ["ede", "hex", "thick", "edethick", "front", "back", "labrmul", "emul", "demul"] .* ".jld" continue end
        println("Saving to $name")
        JLD.save(joinpath(path, name), Dict("value" => if typeof(matrix.matrix) <: Array matrix.matrix else sdata(matrix.matrix) end, "maxx"  => matrix.maxx,
                                            "maxy"  => matrix.maxy,
                                            "title" => matrix.title,
                                            "minx"  => matrix.minx,
                                            "miny"  => matrix.miny,
                                            "xlabel" => matrix.xlabel,
                                            "ylabel" =>matrix.ylabel))
    end
end

function save(matrix::OMatrix{T}, name, path) where T
    name = replace(name, "Δ"=>"d")*".jld"
    JLD.save(joinpath(path, name), Dict("value"=> matrix.matrix,
                                        "maxx"  => matrix.maxx,
                                        "maxy"  => matrix.maxy,
                                        "title" => matrix.title,
                                        "minx"  => matrix.minx,
                                        "miny"  => matrix.miny,
                                        "xlabel" => matrix.xlabel,
                                        "ylabel" =>matrix.ylabel))
end


function calibrate(x, poly)
    return sum(x.^i .* poly[i] for i in 1:length(poly))
end

function calibrate(x, shift, gain, gain²)
    return shift + x*gain + x^2*gain²
end

function calibrate(x::T, shift, gain) where T <: Integer
    return shift + (x+rand()*2-1)*gain
end

function calibrate(x::T, shift, gain) where T <: AbstractFloat
    return shift + (x+rand()*2-1)*gain
end

function ucalibrate(x, shift, gain)
    return shift + x*gain
end

function ucalibrate(x, shift, gain, gain²)
    return shift + x*gain + x^2 * gain²
end

function ucalibrate(x, shift, gain, gain², gain³)
    return shift + x*gain + x^2 * gain² + x^3*gain³
end

function gainshift!(m::Array{T, 1}, shift, gain, gain²) where T <: Integer
    tmp = zeros(T, size(m))
    first = firstindex(m)
    last = lastindex(m)
    @inbounds @fastmath for i in eachindex(m)
        j = unsafe_trunc.(T, calibrate.(i .+ 2 .*rand(m[i]) .- 1, shift, gain, gain²))
        for _j in j[first .≤ j .≤ last]
            tmp[_j] += 1
        end
    end
    m .= tmp
end

function gainshift!(m::Array{T, 1}, shift, gain) where T <: Integer
    tmp = zeros(T, size(m))
    first = firstindex(m)
    last = lastindex(m)
    @inbounds @fastmath for i in eachindex(m)
        j = unsafe_trunc.(T, ucalibrate.(i .+ 2 .*rand(m[i]) .- 1, shift, gain))
        for _j in j[first .≤ j .≤ last]
            tmp[_j] += 1
        end
    end
    m .= tmp
end

function gainshift!(m::OMatrix{Array{T, 1}}, shift, gain, gain²) where T <: Integer
    gainshift!(m.matrix, shift, gain, gain²)
end

function gainshift(m::Array{T, 1}, shift, gain, gain²) where T <: Integer
    tmp = zeros(T, size(m))
    first = firstindex(m)
    last = lastindex(m)
    @inbounds @fastmath for i in eachindex(m)
        j = unsafe_trunc.(T, ucalibrate.(i .+ 2 .*rand(m[i]) .- 1, shift, gain, gain²))
        for _j in j[first .≤ j .≤ last]
            tmp[_j] += 1
        end
    end
    tmp
end

function gainshift(m::Array{T, 1}, shift, gain) where T <: Integer
    tmp = zeros(T, size(m))
    first = firstindex(m)
    last = lastindex(m)
    @inbounds @fastmath for i in eachindex(m)
        j = unsafe_trunc.(T, ucalibrate.(i .+ 2 .*rand(m[i]) .- 1, shift, gain))
        for _j in j[first .≤ j .≤ last]
            tmp[_j] += 1
        end
    end
    tmp
end

function gainshift(m::Array{T, 1}, shift, gain, gain², gain³) where T <: Integer
    tmp = zeros(T, size(m))
    first = firstindex(m)
    last = lastindex(m)
    @inbounds @fastmath for i in eachindex(m)
        j = unsafe_trunc.(T, ucalibrate.(i .+ 2 .*rand(m[i]) .- 1, shift, gain, gain², gain³))
        for _j in j[first .≤ j .≤ last]
            tmp[_j] += 1
        end
    end
    tmp
end

function gainshift(m::OMatrix{Array{T, 2}}, coeffsx, coeffxy) where T <: Integer
    throw("Gainshift for OMatrix not yet implemented")
end

function increment!(matrix::OMatrix{T} where T<:AbstractArray{Int64, 2}, i, j)
    if matrix.minx < i < matrix.maxx && matrix.miny < j < matrix.maxy
        @inbounds matrix.matrix[xtoindex(i, matrix), ytoindex(j, matrix)] += 1
    end
    nothing
end

function increment!(matrix::OMatrix{T} where T<:AbstractArray{Int64, 1}, i)
    if matrix.minx < i < matrix.maxx
        @inbounds matrix.matrix[xtoindex(i, matrix)] += 1
    end
    nothing
end

function xtoindex(x, matrix::OMatrix)::Int64
    return @inbounds unsafe_trunc(Int64, matrix.α[1]*(x-matrix.minx)) + 1
end

function ytoindex(y, matrix::OMatrix)::Int64
    return @inbounds unsafe_trunc(Int64, matrix.α[2]*(y-matrix.miny)) + 1
end

function indextox(i, matrix::OMatrix)::Float64
    return @inbounds i/matrix.α[1] + matrix.minx
end

function indextoy(j, matrix::OMatrix)::Float64
    return @inbounds j/matrix.α[2] + matrix.miny
end

function xrange(matrix::OMatrix)
    range(matrix.minx, matrix.maxx, length=size(matrix.matrix, 1))
end

function yrange(matrix::OMatrix)
    range(matrix.miny, matrix.maxy, length=size(matrix.matrix, 2))
end

function calibrate(x, event::MiniEvent, value::Symbol, parameters::Parameters)::Float64
    if value == :e || value == :Δe
        coeffs = parameters.gainshift[value][event.back, event.front, :]::Array{Float64, 1}
        return calibrate(x, coeffs...)
    elseif value == :ex
        coeffs = parameters.exfromeΔe[event.front, :]::Array{Float64, 1}
        theoretical = calibrate(x, coeffs...) 
        return parameters.excorr[event.front, 1] + theoretical*parameters.excorr[event.front, 2]
    elseif value == :range

    end
    error("Unsupported symbol for calibration")
end

function calibrate(event::MiniEvent, value::Symbol, parameters::Parameters)::Float64
    coeffs = parameters.gainshift[value][event.back, event.front, :]::Array{Float64, 1}
    return calibrate(getfield(event, value), coeffs...)
end
