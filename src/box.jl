import Base: getindex, setindex!, firstindex, lastindex, size

mutable struct Box{T}<:AbstractMatrix{T}
    values::T
    xindices::UnitRange{Int}
    yindices::UnitRange{Int}
    middle::Tuple{Int, Int}
    function Box(matrix::AbstractMatrix, startx::V, stopx::V, starty::V, stopy::V) where V <: Integer
        @assert startx > 0 && stopx < size(matrix, 1)
        @assert starty > 0 && stopy < size(matrix, 2)
        mid = (middle(startx, stopx), middle(starty, stopy))
        xindices = startx:stopx
        yindices = starty:stopy
        new(view(matrix, xindices, yindices), xindices, yindices, middle)
    end
end

function Box(matrix::AbstractMatrix, indices)
    Box(matrix, indices[1]..., indices[2]...)
end

size(box::Box) = size(box.values)
getindex(box::Box, i) = getindex(box.values, i)
setindex!(box::Box, i, v) = setindex!(box.values, i, v)
firstindex(box::Box) = firstindex(box.values)
lastindex(box::Box) = lastindex(box.values)


middle(start, stop) = (stop-start)/2 + start |> x -> round(Int64, x)
function boxlim(start, stop, width)
    lower = middle(start, stop) - width/2
    upper = middle(start, stop) + width/2
    (lower, upper) .|> x -> floor(Int, x)
end

