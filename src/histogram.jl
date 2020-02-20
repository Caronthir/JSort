function binner(edges::AbstractRange)
    start = first(edges)
    stop = last(edges)
    dx = step(edges)
    maxi = length(edges)
    function tobin(x)::Int64
        if x < start
            1
        elseif x > stop
            maxi
        else
            trunc(Int, (x - start)/dx) + 1
        end
    end
    tobin
end

function histogram(X, Y, xedges, yedges)
    tobinx = binner(xedges)
    tobiny = binner(yedges)
    weights = zeros((length(xedges), length(yedges)))
    for (x, y) in zip(X, Y)
        i = tobinx(x)
        j = tobiny(y)
        weights[i, j] += 1
    end
    weights
end

function histogram(X::AbstractVector, edges::AbstractRange)::Vector{Int}
    tobin = binner(edges)
    weights = zeros(length(edges))
    for x in X
        i = tobin(x)
        weights[i] += 1
    end
    weights
end

function histogram_edges(X::AbstractVector; nbins=1000)
    edges = range(minimum(X), stop=maximum(X), length=nbins)
end

function histogram_edges(X, Y; nbins=1000)
    xedges = range(minimum(X), stop=maximum(X), length=nbins)
    yedges = range(minimum(Y), stop=maximum(Y), length=nbins)
    xedges, yedges
end

histogram_edges(X::AbstractMatrix; nbins=1000) = histogram_edges(X[:, 1], X[:, 2], nbins=nbins)
function histogram(X::AbstractMatrix; nbins=1000)
    edges = histogram_edges(X, nbins=nbins)
    histogram(X[:, 1], X[:, 2], edges...), edges
end

function histogram(X::AbstractVector; nbins=1000)
    edges = histogram_edges(X, nbins=nbins)
    histogram(X, edges), edges
end

