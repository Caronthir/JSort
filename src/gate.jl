import Base:∪, ∩, isless, size

abstract type AbstractGate end

struct Gate{T,V} <: AbstractGate
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
Gate(low, high) = Gate(low=low, high=high)
Gate(mean; width=50) = Gate(mean-width, mean+width)


(gate::Gate{T,V} where {T,V})(x)::Bool        = gate.low < x < gate.high
(gate::Gate{Nothing,V} where {V})(x)::Bool    = x < gate.high
(gate::Gate{T,Nothing} where {T})(x)::Bool    = gate.low < x
function (gate::Gate{Nothing,Nothing})(x)::Bool
    return true
    #throw(ArgumentError("Set the gate limit(s)"))
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

center(g::Gate) = g.low + 1/2*(g.high - g.low)
isless(g1::Gate, g2::Gate) = center(g1) < center(g2)

struct Point
    x::Float64
    y::Float64
end
Point(p::Union{AbstractVector, Tuple{Number, Number}}) = Point(p[1], p[2])

struct Gate2D <: AbstractGate
    lowerleft::Point
    lowerright::Point
    higherright::Point
    higherleft::Point
    function Gate2D(p1::Point, p2::Point, p3::Point, p4::Point)
        p = [p1, p2, p3, p4]
        if !isrectangle(p)
            throw(ArgumentError("Points do not constitute a rectangle"))
        end
        lessy(p1, p2) = p1.y < p2.y
        lessynx(p1, p2) = p1.y <= p2.y && p1.x < p2.x
        p = sort(sort(p, lt=lessy), lt=lessynx)
        new(p[1], p[2], p[4], p[3])
    end
end

Gate2D(points::AbstractVector{Point}) = Gate2D(points...)

function Gate2D(ll::T, lr::T, hr::T, hl::T) where T <: Union{AbstractVector, Tuple{Number, Number}}
    Gate2D(Point(ll), Point(lr), Point(hr), Point(hl))
end

function Gate2D(center::Point; width=50)
    d = width/2
    Gate2D((center.x - d, center.y-d),
         (center.x + d, center.y - d),
         (center.x + d, center.y + d),
         (center.x - d, center.y + d))
end

function Gate2D(xmin::Number, xmax::Number, ymin::Number, ymax::Number)
    Gate2D((xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax))
end

(gate::Gate2D)(p::Point) = gate(p.x, p.y)

function (gate::Gate2D)(x, y)::Bool
    if gate.lowerleft.x < x < gate.lowerright.x && gate.lowerleft.y < y < gate.higherleft.y
        true
    else
        false
    end
end

function isrectangle(points::AbstractVector{Point})::Bool
    if length(points) != 4
        return false
    end

    for i in 1:4
        p = points[i]
        found_x = 0
        found_y = 0
        for j in 1:4
            i == j && continue
            p2 = points[j]
            if p2.x ≈ p.x
                found_x += 1
            end
            if p2.y ≈ p.y
                found_y += 1
            end
        end
        if found_x != 1 || found_y != 1
            return false
        end
    end
    true
end

size(gate::Gate2D) = gate.lowerleft.x, gate.lowerright.x, gate.lowerleft.y, gate.higherleft.y
center(gate::Gate2D) = (gate.lowerleft.x + 1/2*(gate.lowerright.x - gate.lowerleft.x),
                        gate.lowerleft.y + 1/2*(gate.higherright.y - gate.lowerleft.y))
