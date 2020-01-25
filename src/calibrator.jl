using .JSort
using Polynomials
import Base: show, read, fill!
import Printf.@sprintf
using DelimitedFiles

abstract type AbstractCalibrator end

struct ParticleCalibrator{T} <: AbstractCalibrator
    poly::Matrix{Poly{T}}
end


function ParticleCalibrator(coefficients=[0.0, 1.0]; var::Symbol=:x, T=nothing)
    if !isnothing(T)
        coefficients = convert(Vector{T}, coefficients)
    end
    poly = [Poly(coefficients, var) for f in 1:8, b in 1:8]
    ParticleCalibrator(poly)
end

ParticleCalibrator{T}(;var::Symbol=:x) where T = ParticleCalibrator(var=var, T=T)

function ParticleCalibrator(coefficients::AbstractMatrix; var::Symbol=:x)
    @assert size(coefficients) == (8, 8)
    poly = [Poly(coefficients[f, b], var) for f in 1:8, b in 1:8]
    ParticleCalibrator(poly)
end

function show(io::IO, calib::ParticleCalibrator)
    for row in 1:size(calib.poly, 1), col in 1:size(calib.poly, 2)
        string = "f:$row b:$col "
        p = calib.poly[row, col]
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
              string *= @sprintf "%.2f" coefficient
          elseif exponent == 1
              string *= @sprintf "%.2f%s" coefficient p.var
          elseif exponent == 2
              string *= @sprintf "%.2f%s²" coefficient p.var
          else
              string *= @sprintf "%.2f%s^%.2g" coefficient p.var exponent
          end

          if exponent+1 < length(p.a)
              string *= " + "
          end
      end
        println(string)
    end
end

function read(path::AbstractString, V::Type{ParticleCalibrator}; var=:x)
    calibrator = open(path, "r") do io
        read(io, V, var=var)
    end
end

function read(io::IO, ::Type{ParticleCalibrator{T}}; var=:x) where T
    body, header = readdlm(io, ',', header=true)
    b = f = 1
    bfs = Set((f, b) for f in 1:8, b in 1:8)
    coefficients = [T[] for f in 1:8, b in 1:8]
    for (i, row) in enumerate(eachrow(body))
        f = row[1] |> Int
        b = row[2] |> Int
        pop!(bfs, (f, b))
        coefficients[f, b] = row[3, :]
    end
    if length(bfs) > 0
        @warn "Coefficients file lacked " bfs
    end
    ParticleCalibrator(coefficients, var=var)
end

read(io::IO, ::Type{ParticleCalibrator}; var=var) = read(io, ParticleCalibrator{Float32}, var=var)

function (calib::ParticleCalibrator)(front::Integer, back::Integer, x)
    calib.poly[front, back](x)
end

function fill!(calib::ParticleCalibrator, x)
    fill!(calib.poly, x)
end

function combine!(destination::ParticleCalibrator, source::ParticleCalibrator)
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
