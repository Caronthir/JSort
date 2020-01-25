

function calibrate(x, poly)
    return sum(x.^i .* poly[i] for i in 1:length(poly))
end

function calibrate(x, shift, gain, gain²)
    return shift + x*gain + x^2*gain²
end

function calibrate(x::Real, shift, gain)
    return shift + (x+rand()*2-1)*gain
end

function ucalibrate(x, shift, gain)
    return shift + x*gain
end

function ucalibrate(x, shift, gain, gain²)
    return shift + gain*x + gain²*x^2
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
