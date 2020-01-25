using .JSort
import Base:read
import Base:show
using DataFrames

struct ParticleDetector{T, V}
    e::Vector{T}
    Δe::Vector{V}
    front::Int8
    back::Int8
end

function ParticleDetector(e::AbstractArray{T}, Δe::AbstractArray{V}, f::Number, b::Number) where {T, V}
    ParticleDetector(e |> vec, Δe |> vec, f, b)
end

readparticledetector(path::AbstractString; T=Float32, V=Float32) = read(path, ParticleDetector{T, V})

function read(path::AbstractString, V::Type{ParticleDetector})
    detector = open(path, "r") do io
        read(io, V)
    end
end

function read(io::IO, ::Type{ParticleDetector{T, T}}) where T
    # Assumes binary format
    len = read(io, T) |> Int
    data = Matrix{T}(undef, len, 2)
    read!(io, data)
    ParticleDetector(data[:, 1], data[:, 2])
end

read(io::IO, ::Type{ParticleDetector}) = read(io, ParticleDetector{Float32, Float32})

function show(io::IO, detector::ParticleDetector)
    println("Front: $(detector.front)\nBack: $(detector.back)")
    show(DataFrame(e=detector.e, Δe=detector.Δe), allrows=false)
end

function (calib::ParticleCalibrator)(d::ParticleDetector)
    poly = calib[d.front, d.back]
    if poly.var == :e
        map!(poly, d.e, d.e)
    elseif poly.var == :Δe
        map!(poly, d.Δe, d.Δe)
    else
        throw("Unable to deduce how to apply calibration. Specify `e` or `Δe`")
    end
end

function (calib::ParticleCalibrator)(d::ParticleDetector, var::Symbol)
    poly = calib[d.front, d.back]
    if var == :e
        map!(poly, d.e, d.e)
    elseif var == :Δe
        map!(poly, d.Δe, d.Δe)
    else
        throw("Variable must be `e` or `Δe`, not $var")
    end
end

# function addcalibrator!(detector::ParticleDetector, calib::ParticleCalibrator, ::Type{Val{:e}})
#     combine!(detector.ecalibrator, calib)
# end

# function addcalibrator!(detector::ParticleDetector, calib::ParticleCalibrator, ::Type{Val{:Δe}})
#     combine!(detector.Δecalibrator, calib)
# end

# function addcalibrator!(detector::ParticleDetector, calib::ParticleCalibrator)
#     addcalibrator!(detector::ParticleDetector, calib::ParticleCalibrator, calib.poly[1].var)
# end

# addcalibrator!(detector::ParticleDetector, calib::ParticleCalibrator, s::Symbol) = addcalibrator!(detector, calib, Val{s})

# function calibrate!(pd::ParticleDetector)
#     #=
#       YOU MUMBLING IDIOT! WRONG TYPE STRUCTURE!!! AAARGH!

#     =#
#     for i in eachindex(pd.e)
#         pd.e[i] = pd.ecalibrator(pd.e[i])
#     end
#     for i in eachindex(pd.Δe)
#         pd.Δe[i] = pd.Δecalibrator(pd.Δe)
#     end
# end
