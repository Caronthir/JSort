using .JSort
import Base: isless

struct GammaCascade
    gate::AbstractGate
    e::Vector{Float32}
    regions::Vector{Tuple{Float32, Float32}}
    keywords::Dict{Symbol, Any}
    #function GammaCascade(gate::Gate, e::AbstractVector, regions::AbstractVector)
        #new(gate, sort(e) |> reverse)
    #end
end
GammaCascade(gate::AbstractGate, e::AbstractVector, regions::AbstractVector) = GammaCascade(gate, e, regions, Dict(:align => Dict()))

function extend!(base::GammaCascade, others::AbstractArray{GammaCascade})
    lower = [c for c in sort(others) if c < base]
    for other in lower
        extend!(base, other)
    end
end

function extend!(base::GammaCascade, other::GammaCascade; width=3)
    from = center(base.gate)
    to = center(other.gate)
    for e in base.e
        if abs(from - e - to) < width
            append!(base.e, other.e)
            break
        end
    end
end

isless(c1::GammaCascade, c2::GammaCascade) = isless(c1.gate, c2.gate)
center(c::GammaCascade) = center(c.gate)
