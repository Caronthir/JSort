using JSort
using PyCall
using Statistics

function featurealign2d(reference::OMatrix{T}, target::OMatrix{T}; numregions=5, width=nothing,
                        searchwidth=50, smoother=nothing, limit=50) where T
    featurealign2d(reference.matrix, target.matrix, numregions=numregions, width=width,
                   searchwidth=searchwidth, smoother=smoother, limit=limit)
end

function featurealign2d(reference, target; numregions=5, width=nothing,
                        searchwidth=50, smoother=nothing, limit=50)
    if !isnothing(smoother)
        throw("Smoothing not implemented")
    end
    subarrays = splitregions2d(reference, numregions, width=width, limit=limit)

    targetpeaks_x = Int64[]
    targetpeaks_y = Int64[]
    referencepeaks_x = Int64[]
    referencepeaks_y = Int64[]
    for array in subarrays
        tx, ty = featureoverlap2d(array, target, searchwidth)
        rx, ry = array.indices
        rm = (rx, ry) .|> middle
        tm = (tx, ty) .|> middle
        push!(targetpeaks_x, tm[1])
        push!(targetpeaks_y, tm[2])
        push!(referencepeaks_x, rm[1])
        push!(referencepeaks_y, rm[2])
    end
    referencepeaks_x, targetpeaks_x, referencepeaks_y, targetpeaks_y
end


function featureoverlap2d(feature, target, searchwidth)
    x, y = feature.indices 
    widthx, widthy = size(feature)
    startx  = clip(first(x) - searchwidth, size(target, 1))
    stopx   = clip(last(x)  + searchwidth, size(target, 1))
    starty  = clip(first(y) - searchwidth, size(target, 2))
    stopy   = clip(last(y)  + searchwidth, size(target, 2))
    windowx = startx:(startx + widthx - 1)
    windowy = starty:(starty + widthy - 1)
    lengthx, lengthy = (stopx - startx - widthx,
                        stopy - starty - widthy)
    overlap = zeros((lengthx, lengthy))
    for i in 1:lengthx, j in 1:lengthy
        shiftx = i - 1  # Account for zero shift
        shifty = j - 1
        diff = (feature .- target[windowx .+ shiftx,
                                  windowy .+ shifty]).^2
        overlap[i, j] = sum(diff)
    end
    return overlap

    # Find the best overlap
    indices = argmin(overlap)
    x = range(startx + indices[1] + 1, length=widthx)
    y = range(starty + indices[2] + 1, length=widthy)
    (x, y)
end


function arraybox(array, startx::T, stopx::T, starty::T, stopy::T) where T <: Integer
    subarray = view(array, startx:stopx, starty:stopy)
end

function splitregions2d(array::OMatrix{T}, numberofregions; width=nothing, limit=50) where T
    splitregions2d(array, numberofregions, width=width, limit=limit)
end

function splitregions2d(array, numberofregions; width=nothing, limit=50)
    edges = boxbounds(array, numberofregions, width=width)
    arrays = []
    for edge in edges
        subarray = arraybox(array, edge)
        if sum(subarray) > limit
            push!(arrays, subarray) 
        end
    end
    arrays
end

function arraybox(array, edges)
    arraybox(array, edges[1]..., edges[2]...)
end

function boxbounds(array, numberofregions; width=nothing)
    sizex, sizey = size(array)
    stepx, stepy = floor.(Int, size(array)./numberofregions)
    if isnothing(width)
        width = stepx
    end
    edges = []
    for x in 1:numberofregions
        startx = (x-1)stepx + 1
        stopx = x*stepx + 1
        lowerx, upperx = boxlim(startx, stopx, width) .|> x -> clip(x, sizex)
        for y in 1:numberofregions
            starty = (y-1)stepy + 1
            stopy = y*stepy + 1
            lowery, uppery = boxlim(starty, stopy, width) .|> y -> clip(y, sizey)
            push!(edges, ((lowerx, upperx), (lowery, uppery))) 
        end
    end
    edges
end

middle(start, stop) = (stop-start)/2 + start |> x -> round(Int64, x)
function boxlim(start, stop, width)
    lower = middle(start, stop) - width/2
    upper = middle(start, stop) + width/2
    (lower, upper) .|> x -> floor(Int, x)
end
