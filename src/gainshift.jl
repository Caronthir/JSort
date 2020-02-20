

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
