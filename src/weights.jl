abstract type Weight end

struct ThresholdSum <: Weight
    threshold::Float64
    function ThresholdSum(threshold::Real)
        0 < threshold && return new(1.0 + threshold)
        throw(ArgumentError("Threshold must be less than 1"))
    end
end

struct ThresholdArea <: Weight
    threshold::Float64
    function ThresholdArea(threshold::Real)
        0 < threshold && return new(1.0 + threshold)
        throw(ArgumentError("Threshold must be less than 1"))
    end
end

function (w::ThresholdSum)(x)
    min = minimum(x)
    sum(x .≤ w.threshold*min)
end

function (w::ThresholdArea)(x)
    min = minimum(x)
    thresholded = x .≤ w.threshold*min
    maxx = maxy = -1
    miny = minx = 1e5
    for x in 1:size(thresholded, 1), y in 1:size(thresholded, 2)
        if thresholded[x, y]
            if x < minx
                minx = x
            end
            if x > maxx
                maxx = x
            end
            if y < miny
                miny = y
            end
            if y > maxy
                maxy = y
            end
        end
    end
    sidex = maxx - minx
    sidey = maxy - miny
    (sidex > 0 ? sidex : 1) * (sidey > 0 ? sidey : 1)
end
