using .JSort
using Polynomials
import Base: show, read, fill!, write, getindex, setindex!
import Printf.@sprintf
using DelimitedFiles: readdlm, writedlm

abstract type AbstractCalibrator end
abstract type AbstractCalibrator2D <: AbstractCalibrator end
abstract type AbstractCalibrator1D <: AbstractCalibrator end

function show(io::IO, calib::AbstractCalibrator1D)
    string = ""
    for (chn, p) in enumerate(calib.poly)
        string *= "Channel $chn:\t" * serializepoly(p)
        string *= "\n"
    end
    print(io, string)
end

function show(io::IO, calib::AbstractCalibrator2D)
    string = ""
    for row in 1:size(calib.poly, 1), col in 1:size(calib.poly, 2)
        string *= "f:$row b:$col\t"
        p = calib.poly[row, col]
        string *= serializepoly(p) * "\n"
    end
    print(io, string)
end

function read(path::AbstractString, ::Type{T}; var=:x) where T<:AbstractCalibrator
    calibrator = open(path, "r") do io
        read(io, T, var=var)
    end
end

function read(io::IO, ::Type{T}; var=:x) where T<:AbstractCalibrator1D
    body, header = readdlm(io, ',', header=true)
    channel = 1
    chns = Set(chn for chn in 1:33)
    polys = Poly[Poly{Float32}([0.0]) for i in 1:33]
    for (i, row) in enumerate(eachrow(body))
        chn = row[1] |> Int
        pop!(chns, chn)
        coeff = filter(x -> length(x) > 0, row[2:end]) .|> Float32
        polys[chn] = Poly(coeff, var)
    end
    if length(chns) > 0
        @warn "Gamma coefficients file lacked channel:\n" chns
    end
    T(polys)
end

function read(io::IO, ::Type{T}; var=:x) where T<:AbstractCalibrator2D
    body, header = readdlm(io, ',', header=true)
    b = f = 1
    bfs = Set((f, b) for f in 1:8, b in 1:8)
    #coefficients = [Float32[] for f in 1:8, b in 1:8]
    polys = Poly[Poly{Float32}([0.0]) for f in 1:8, b in 1:8]
    for (i, row) in enumerate(eachrow(body))
        f = row[1] |> Int
        b = row[2] |> Int
        pop!(bfs, (f, b))
        #coefficients[f, b] = row[3:end]
        coeff = filter(x -> length(x) > 0, row[3:end]) .|> Float32
        polys[f, b] = Poly(coeff, var)
    end
    if length(bfs) > 0
        @warn "Coefficients file lacked (f, b):\n" bfs
    end
    T(polys)
end

function write(path::AbstractString, calib::AbstractCalibrator1D)
    lines = Any[["channel", ["x^$(i-1)" for i in 1:length(calib.poly[1, 1].a)]...]]
    for (chn, poly) in enumerate(calib.poly)
        push!(lines, [chn, poly.a...])
    end
    writedlm(path, lines, ',')
end

function write(path::AbstractString, calib::AbstractCalibrator2D)
    lines = Any[["f", "b", ["x^$(i-1)" for i in 1:length(calib.poly[1, 1].a)]...]]
    for f in 1:size(calib.poly, 1), b in 1:size(calib.poly, 2)
        poly = calib.poly[f, b]
        push!(lines, [f, b, poly.a...])
    end
    writedlm(path, lines, ',')
end

function serializepoly(p::Poly)::String
    string = ""
    for (exponent, coefficient) in enumerate(p.a)
        exponent -= 1

        # Special case of p = 0.0
        if coefficient == 0
            if exponent == 0 && length(p.a) == 1
                string *= "0"
            end
            continue
        end

        if exponent == 0
            string *= @sprintf "%6.2e" coefficient
        elseif exponent == 1
            string *= @sprintf "%5.2e%s" coefficient p.var
        elseif exponent == 2
            string *= @sprintf "%.2e%s²" coefficient p.var
        else
            string *= @sprintf "%.2e%s^%.2g" coefficient p.var exponent
        end

        if exponent+1 < length(p.a)
                string *= " +\t"
        end
    end
    return string
end

fill!(calib::AbstractCalibrator, x) = fill!(calib.poly, x)

function combine!(destination::AbstractCalibrator, source::AbstractCalibrator)
    @assert size(destination.poly) == size(source.poly)
    @assert destination.poly[1].var == source.poly[1].var
    var = destination.poly[1].var
    for i in eachindex(destination.poly)
        # Only works for var = :x ???
        p = Poly(source.poly[i].a, :x)
        q = Poly(destination.poly[i].a, :x)
        destination.poly[i] = Poly(q(p).a, var)
    end
end

invert(calib::AbstractCalibrator) = typeof(calib)(map(invert, calib.poly))
function invert(poly::Poly)
    degree = length(poly.a) - 1
   if degree == 0
        Poly([-poly.a], poly.var)
    elseif degree == 1
        b, a = poly.a
        Poly([-b/a, 1/a], poly.var)
    else
        error("Polynomial inversion only supported for degree 0 and 1.")
    end
end

getindex(X::AbstractCalibrator, i) = getindex(X.poly, i)
getindex(X::AbstractCalibrator, I...) = getindex(X.poly, I...)
setindex!(X::AbstractCalibrator, i) = setindex!(X.poly, i)
setindex!(X::AbstractCalibrator, I...) = setindex!(X.poly, I...)

(calib::AbstractCalibrator)(x::Number) = map(y -> y(x), calib.poly)

## Particle Calibrator
struct ParticleCalibrator <: AbstractCalibrator2D
    poly::Matrix{Poly}
end

function ParticleCalibrator(coefficients::AbstractVector{<:Number}=[0.0, 1.0]; var::Symbol=:x, T=nothing)
    if !isnothing(T)
        coefficients = convert(Vector{T}, coefficients)
    end
    poly = [Poly{Float32}(coefficients, var) for f in 1:8, b in 1:8]
    ParticleCalibrator(poly)
end

ParticleCalibrator(path::AbstractString; var::Symbol=:x) = read(path, ParticleCalibrator, var=var)

function ParticleCalibrator(coefficients::AbstractMatrix{<:Number}; var::Symbol=:x)
    @assert size(coefficients) == (8, 8)
    poly = [Poly{Float32}(coefficients[f, b], var) for f in 1:8, b in 1:8]
    ParticleCalibrator(poly)
end


#read(io::IO, ::Type{ParticleCalibrator}; var=:x) = read(io, ParticleCalibrator, var=var)

## Time calibration
struct TimeCalibrator <: AbstractCalibrator1D
    poly::Vector{Poly{Float32}}
end
#TimeCalibrator(p::Vector{<:Poly}) = TimeCalibrator()

function TimeCalibrator(coefficients::AbstractMatrix{<:Number})
    poly = [Poly(coeff, :t) for coeff in eachrow(coefficients)]
    TimeCalibrator(poly)
end

function TimeCalibrator(coefficients::AbstractVector{<:Number}=[0.0, 1.0]; var::Symbol=:t)
    poly = [Poly(coefficients, :t) for chn in 1:33]
    TimeCalibrator(poly)
end
TimeCalibrator(path::AbstractString) = read(path, TimeCalibrator, var=:t)

## Time Correction
struct TimeFunction
    t₀::Float32
    c1::Float32
    c2::Float32
    c3::Float32
end
TimeFunction() = TimeFunction(0,0,0,0)
function (tf::TimeFunction)(E)
    return tf.t₀ - tf.c1/(E - tf.c2) + tf.c3*E
end
function (tf::TimeFunction)(t, E)
    return t - tf.t₀ - tf.c1/(E - tf.c2) + tf.c3*E
end

function show(io::IO, tf::TimeFunction)
    s = @sprintf "%.2e + %.2e/(E + %.2e)⁻¹ + %.2e×E" tf.t₀ tf.c1 tf.c2 tf.c3
    println(io, s)
end

struct TimeCorrector <: AbstractCalibrator1D
    poly::Vector{TimeFunction}
end
TimeCorrector(;var::Symbol=:t) = TimeCorrector([TimeFunction() for channel in 1:33])
TimeCorrector(path::AbstractString) = read(path, TimeCorrector)

function show(io::IO, tc::TimeCorrector)
    text = ""
    for (chn, tf) in enumerate(tc.poly)
        text *= "Channel $chn:\t"
        text *= repr(tf)
    end
    println(text)
end

function read(io::IO, ::Type{TimeCorrector}; var=:t)
    body, header = readdlm(io, ',', header=true)
    chn = 1
    chns = Set(channel for channel in  1:33)
    time = TimeFunction[TimeFunction()]
    for (i, row) in enumerate(eachrow(body))
        @assert length(row) == 5 "Corrupt time at row $i"
        chn = row[1] |> Int
        pop!(chns, chn)
        time[chn] = TimeFunction(row[2:end]...)
    end
    if length(chns) > 0
        @warn "Time file lacked channel(s):\n" chns
    end
    TimeCorrector(time)
end

function write(path::AbstractString, calib::TimeCorrector)
    lines = Any[["channel", "t₀", "c₁", "c₂", "c₃"]]
    for (chn, tf) in enumerate(calib.poly)
        push!(lines, [chn, tf.t₀, tf.c1, tf.c2, tf.c3])
    end
    writedlm(path, lines, ',')
end

## Gamma calibrator
struct GammaCalibrator <: AbstractCalibrator1D
    poly::Vector{Poly{Float32}}
end
function GammaCalibrator(coefficients::AbstractMatrix{<:Number})
    poly = [Poly(coeff, :eᵧ) for coeff in eachrow(coefficients)]
    GammaCalibrator(poly)
end

function GammaCalibrator(coefficients::AbstractVector{<:Number}=[0.0, 1.0]; var::Symbol=:eᵧ)
    poly = [Poly(coefficients, :eᵧ) for chn in 1:33]
    GammaCalibrator(poly)
end
GammaCalibrator(path::AbstractString) = read(path, GammaCalibrator, var=:eᵧ)
