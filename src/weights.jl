abstract type Weight end

#=
    Sum over the bins that are within `threshold`% of the minimum.
    A large value means the well is shallow.
=#
struct ThresholdSum <: Weight
    threshold::Float64
    function ThresholdSum(threshold::Real)
        0 < threshold && return new(1.0 + threshold)
        throw(ArgumentError("Threshold must be less than 1"))
    end
end

#=
  Sum over the bins that are within the rectangle inscribing the
  the bins which are `threshold`% of the minimum.
  A large value means the well is shallow. Punished asymmetric wells
  more than ThresholdSum
=#
struct ThresholdArea <: Weight
    threshold::Float64
    function ThresholdArea(threshold::Real)
        0 < threshold && return new(1.0 + threshold)
        throw(ArgumentError("Threshold must be less than 1"))
    end
end

struct ThresholdSumShift <: Weight
    threshold::Float64
    shift::Float64
    function ThresholdSumShift(threshold::Real)
        0 < threshold && return new(1.0 + threshold, 1.0)
        throw(ArgumentError("Threshold must be less than 1"))
    end
end

#=
   Sum over the bins within the `radius` distance of the
   minimum, minus the minimum. A high score means the well is sharp.
=#
struct WellSum <: Weight
    radius::Int
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

function (w::ThresholdSumShift)(x)
    min = minimum(x) + w.shift
    sum(x .≤ w.threshold*min)
end

function (w::WellSum)(x)
    truth = argmin(x)
    r = w.radius
    xstart, ystart = truth.I .- r
    xstop, ystop = truth.I .+ r
    xstart = xstart > 0 ? xstart : 1
    ystart = ystart > 0 ? ystart : 1
    xstop = xstop < size(x, 1) ? xstop : size(x, 1)
    ystop = ystop < size(x, 2) ? ystop : size(x, 2)
    sum(@view x[xstart:xstop, ystart:ystop]) - minimum(x)
end
