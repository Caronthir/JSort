function L1(x, y)
    return sum(abs.(x .- y))
end

function L2(x, y)
    return sum((x .- y).^2)
end

function χ²(x, y)
    S = Float64.((x .- y).^2)
    S[x .!= 0] .= S[x .!= 0]./x[x .!= 0]
    return sum(S)
end

function L1l(x, y, limit)
    sum(abs.((x .≥ limit) .- (y .≥ limit)))
end

function L2l(x, y, limit)
    sum(((x .≥ limit) .- (y .≥ limit)).^2)
end

function splitmetric(x, y; metric=L1, threshold=0.5)
    projection = sum(x, dims=2) |> vec
    sp = argthreshold(projection, threshold)
    # The regions should weigh `threshold` and `(1-threshold)`
    # Add penalization equal to their difference
    l1 = metric(x[:, 1:sp],   y[:, 1:sp])
    l2 = metric(x[:, sp:end], y[:, sp:end])
    diff = abs(threshold*l1 - (1-threshold)*l2)
    return l1 + l2 + diff
end

#=
Return the index where the cumulative sum of `x` is
`threshold` of the total sum.
=#
function argthreshold(x, threshold=0.5)
    csum = cumsum(x)
    limit = threshold*csum[end]
    for (i, x_) in enumerate(csum)
        if x_ ≥ limit
            return i
        end
    end
end
