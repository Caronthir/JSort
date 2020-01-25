import Base:∪, ∩

struct Gate{T,V}
    low::T
    high::V
    function Gate(;low::T=nothing, high::V=nothing) where {T,V}
        if !isnothing(low) && !isnothing(high)
            @assert low < high "Lower limit must be lower than higher: $low ≮ $high"
        end
        new{T,V}(low, high)
    end
end
#Gate(;low=nothing, high=nothing) = Gate(low, high)


(gate::Gate{T,V} where {T,V})(x)::Bool        = gate.low < x < gate.high
(gate::Gate{Nothing,V} where {V})(x)::Bool    = x < gate.high
(gate::Gate{T,Nothing} where {T})(x)::Bool    = gate.low < x
function (gate::Gate{Nothing,Nothing})(x)::Nothing
    throw(ArgumentError("Set the gate limit(s)"))
end

∩(g1::Gate{Nothing, Nothing}, g2::Gate{Nothing, Nothing}) = g1
function ∩(g1::Gate, g2::Gate)
    low = if !isnothing(g1.low) && !isnothing(g2.low)
        g1.low > g2.low ? g1.low : g2.low
        elseif !isnothing(g1.low)
            g1.low
        else
            g2.low
    end
    high = if !isnothing(g1.high) && !isnothing(g2.high)
        g1.high < g2.high ? g1.high : g2.high
        elseif !isnothing(g1.high)
            g1.high
        else
            g2.high
    end
    Gate(low=low, high=high)
end


