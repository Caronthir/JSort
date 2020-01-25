using JSort
using Statistics
using Distributions: Poisson, rand
import PyPlot; const plt = PyPlot
import PyCall
using ProgressBars
using ImageFiltering: imfilter, imfilter!, Kernel
using CSV: read
using DataFrames
using Polynomials
using Base.Filesystem: splitext

function featurealign2d(path::String, pattern::Regex=r"edeb(?<b>\d)f(?<f>\d).bin", referenceid=1
                        ;weight::Weight=ThresholdSum(.3), metric::Function=L2, savepath::String="../figures/",
                        neighbourwidth=9, neighbourthreshold=80, plot=false)
    # Sort all back detectors corresponding to common front detectors
    frontdetectors = Dict(front => Dict{Int,String}() for front in 1:8)
    # Find all files
    anymatch = false
    for fname in readdir(path)
        m = match(pattern, fname)
        if m !== nothing
            f, b = parse.(Int, (m[:f], m[:b]))
            frontdetectors[f][b] = fname
            anymatch = true
        end
    end

    if !anymatch
        throw("Found no matching files using $pattern in $path")
    end

    model = [(x, y, w) -> leastsquares(x, y;order= 1)]
    name = ["Linear Regression"]

    lm = pyimport("sklearn.linear_model")
    loss = ["absolute_loss", "squared_loss"]
    models = [(x, y, w) -> begin
                  X = reshape(x, :, 1)
                  lr.fit(X, y, w)
                  return [lr.estimator_.intercept_; lr.estimator_.coef_]
              end
              for lr in [lm.RANSACRegressor(loss=l) for l in loss]]
    names = ["RANSAC "*l for l in loss]
	models = [model..., models...]
	names = [name..., names...]

    coefficientsx = []
    coefficientsy = []

    function readfile(fname::AbstractString)
        path_ = joinpath(path, fname)
        if splitext(fname)[2] == ".bin"
            readede(path_) |> x -> DataFrame(x, [:e, :Δe])
        elseif splitext(fname)[2] == ".csv"
            CSV.read(path_, header=["e", "Δe"], datarow=2) |> DataFrame
        end
    end


    for (front, backdetectors) in frontdetectors
        reference_fname = backdetectors[referenceid]
        reference_df = readfile(reference_fname)
        edges = histogram_edges(reference_df[!, :e], reference_df[!, :Δe], nbins=1000)
        reference_raw = JSort.histogram(reference_df[!, :e], reference_df[!, :Δe], edges...)
        fig, ax = plt.subplots()
        ax.matshow(reference_raw)
        reference = neighbourreduce(reference_raw, width=neighbourwidth, threshold=neighbourthreshold)
        fig, ax = plt.subplots()
        ax.matshow(reference)
        total = copy(reference_raw)
        for (back, fname) in backdetectors
            if back == referenceid
                push!(coefficientsx, [front, back, 0.0, 1.0])
                push!(coefficientsy, [front, back, 0.0, 1.0])
                continue
            end
            target_df = readfile(fname)
            target_raw = JSort.histogram(target_df[!, :e], target_df[!, :Δe], edges...)
            target = neighbourreduce(target_raw, width=neighbourwidth, threshold=80)

            points = featurealign2d_weighted(reference, target, numregions=15,
                                             limit=1000, searchwidth=40, weight=weight,
                                             metric=metric)
            @show points
            rpx, tpx, rpy, tpy, weights = points
            annotate_plot_2d(target, tpx, tpy, weights)
            plt.show()
            return
            # Convert from index basis to energy basis
            rex, tex = edges[1][rpx], edges[1][tpx]
            rey, tey = edges[2][rpy], edges[2][tpy]

			optimal = ()
			for (name, model) in zip(names, models)

                #coeffx, Σx = leastsquares(tex, rex, weights; order=1)
                #coeffy, Σy = leastsquares(tey, rey, weights; order=1);
				coeffx = model(tex, rex, weights)
				coeffy = model(tey, rey, weights)

				px = Poly(coeffx)
				py = Poly(coeffy)

                e =  px.(target_df[!, :e])
                de = py.(target_df[!, :Δe])

                calibrated = JSort.histogram(e, de, edges...);

                err_un = compare_2d(reference, target, metric=metric)
                err_ca = compare_2d(reference, calibrated, metric=metric)
				ratio = err_ca/err_un
				println("Error of $name: $ratio")

				if length(optimal) == 0
                    optimal =  (name, err_ca, err_un, coeffx, coeffy, calibrated)
				elseif optimal[2] > err_ca
					optimal = (name, err_ca, err_un, coeffx, coeffy, calibrated)
                end
			end
            println("Optimal is $(optimal[1]) with $(optimal[2])")
            name, err_ca, err_un, coeffx, coeffy, calibrated = optimal
            push!(coefficientsx, [front, back, coeffx...])
            push!(coefficientsy, [front, back, coeffy...])
            total .+= calibrated
            if plot
              fig, ax = plt.subplots(ncols=2, nrows=2, sharex=true, sharey=true,
                                    figsize=(20, 20))
              # ax[1].matshow(transpose(reference) .|> log10)
              # ax[1].scatter(rpx, rpy, marker="x", c="r", s=0.5)
              # ax[1].invert_yaxis()
              annotate_plot_2d(reference_raw, rpx, rpy, ax=ax[1])
              ax[1].set_title("f$front b$referenceid")
              # ax[2].matshow(transpose(target) .|> log10)
              # ax[2].scatter(tpx, tpy, marker="x", c="r", s=0.5)
              annotate_plot_2d(target_raw, tpx, tpy, weights, ax=ax[2])
              ax[2].set_title("f$front b$back")
              ax[2].invert_yaxis()
              ax[3].matshow(transpose(reference_raw) .|> log10, cmap="Reds")
              ax[3].matshow(transpose(target_raw) .|> log10, cmap="Blues")
              ax[3].annotate(string(err_un), xy=(1,0), xycoords="axes fraction",
                            horizontalalignment="right", verticalalignment="bottom")
              ax[3].set_title("Original")
              ax[3].invert_yaxis()
              ax[4].matshow(transpose(reference_raw) .|> log10, cmap="Reds")
              ax[4].matshow(transpose(calibrated) .|> log10, cmap="Blues")
              ax[4].set_title("Calibrated - "*name)
              ax[4].invert_yaxis()
              ax[4].annotate(string(err_ca), xy=(1,0), xycoords="axes fraction",
                            horizontalalignment="right", verticalalignment="bottom")
              fig.savefig(joinpath(savepath, "medef$(front)b$(back).png"), bbox_inches="tight", dpi=196)
              plt.close();

              fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))
              plotresiduals(rpx, tpx, coeffx, ax=ax[1:2])
              plotresiduals(rpy, tpy, coeffy, ax=ax[3:4])
              fig.savefig(joinpath(savepath, "residuef$(front)b$(back).png"), bbox_inches="tight", dpi=196)
              plt.close();
            end
        end
        if plot
          fig, ax = plt.subplots()
          ax.matshow(total' .|> log10)
          ax.invert_yaxis()
          fig.savefig(joinpath(savepath, "medef$front.png"), bbox_inches="tight", dpi=196)
          plt.close();
        end
    end
    coefficientsx, coefficientsy
end

function featurealign2d(reference, target; numregions=5, width=nothing,
                        searchwidth=50, smoother=nothing, limit=50)
    if !isnothing(smoother)
        throw("Smoothing not implemented")
    end
    reference = copy(reference)
    target = copy(target)
    subarrays = splitregions2d(reference, numregions, width=width, limit=limit)
    imfilter!(target, target, Kernel.gaussian(3))
    target .= target ./ sum(target)
    imfilter!(reference, reference, Kernel.gaussian(3))
    reference .= reference ./ sum(reference)

    targetpoints_x    = Int64[]
    targetpoints_y    = Int64[]
    referencepoints_x = Int64[]
    referencepoints_y = Int64[]
    i = 0
    for array in subarrays
        x, y = featureoverlap2d(array, target, searchwidth)
        rx, ry = array.indices
        rm = (rx, ry) .|> middle

        push!(targetpoints_x, x)
        push!(targetpoints_y, y)
        push!(referencepoints_x, rm[1])
        push!(referencepoints_y, rm[2])
    end
    referencepoints_x, targetpoints_x, referencepoints_y, targetpoints_y
end

function featurealign2d_weighted(reference, target; numregions=5, width=nothing,
                                 searchwidth=50, smoother=nothing, limit=50, 
                                 metric::Function=L2,
                                 weight::Weight=ThresholdArea(.2))
    if !isnothing(smoother)
        throw("Smoothing not implemented")
    end
    reference = copy(reference)
    target = copy(target)
    subarrays = splitregions2d(reference, numregions, width=width, limit=limit)
    imfilter!(target, target, Kernel.gaussian(3))
    target .= target ./ sum(target)
    imfilter!(reference, reference, Kernel.gaussian(3))
    reference .= reference ./ sum(reference)

    targetpoints_x    = Int64[]
    targetpoints_y    = Int64[]
    referencepoints_x = Int64[]
    referencepoints_y = Int64[]
    weights           = Float64[]
    i = 0
    for array in subarrays
        x, y, w = featureoverlap2d_weighted(array, target, searchwidth, weight, metric)
        rx, ry = array.indices
        rm = (rx, ry) .|> middle

        push!(targetpoints_x, x)
        push!(targetpoints_y, y)
        push!(referencepoints_x, rm[1])
        push!(referencepoints_y, rm[2])
        push!(weights, w)
    end
    referencepoints_x, targetpoints_x, referencepoints_y, targetpoints_y, weights
end


function featurealign2d_ensemble(reference, target; numregions=5, width=nothing,
                        searchwidth=50, smoother=nothing, limit=50)
    if !isnothing(smoother)
        throw("Smoothing not implemented")
    end
    reference = copy(reference)
    target = copy(target)
    subarrays = splitregions2d(reference, numregions, width=width, limit=limit)
    target = imfilter(target, Kernel.gaussian(3))
    target .= target ./ sum(target)
    N = sum(reference)

    targetpoints_x    = Int64[]
    targetpoints_y    = Int64[]
    referencepoints_x = Int64[]
    referencepoints_y = Int64[]
    i = 0
    for array in subarrays
        #tx, ty = featureoverlap2d(array, target, searchwidth)
        featureoverlap2d_ensemble(array, target, searchwidth, N)
        # println(ty, tx)
        # rx, ry = array.indices
        # rm = (rx, ry) .|> middle
        # weight by difference in region total sum?
        if (i += 1) > 1
            return
        end

        # push!(targetpoints_x, tm[1])
        # push!(targetpoints_y, tm[2])
        # push!(referencepoints_x, rm[1])
        # push!(referencepoints_y, rm[2])
    end
    referencepoints_x, targetpoints_x, referencepoints_y, targetpoints_y
end

function featureoverlap2d(feature, target, searchwidth, metric::Function=L2)
    window, _ = shiftsinrange(feature, target, searchwidth)
    distance = slidingwindow(feature, target, metric, searchwidth)
    min = argmin(distance)
    fig, ax = plt.subplots()
    ax.matshow(distance .|> log10)
    ax.contour(distance)
    ax.scatter(min[2]-1, min[1]-1, marker="x", c="r")
    # fig, ax = plt.subplots(ncols=2)
    # ax[1].matshow(feature)
    # X, Y = window[1] .+ min[1], window[2] .+ min[2]
    # ax[2].matshow(target[X, Y])
    # Convert from surface coordinates to feature indices
    offset = first.(window) .+ (length.(window) .÷ 2)
    return min[1] + offset[1], min[2] + offset[2]
end

function featureoverlap2d_weighted(feature, target, searchwidth, weight::Weight, metric::Function=L2)
    window, _ = shiftsinrange(feature, target, searchwidth)
    distance = slidingwindow(feature, target, metric, searchwidth)
    min = argmin(distance)
    # fig, ax = plt.subplots()
    # ax.matshow(distance .|> log10)
    # ax.contour(distance)
    # ax.scatter(min[2]-1, min[1]-1, marker="x", c="r")
    # fig, ax = plt.subplots()
    # d = distance .<= 1.2.*distance[min]
    # ax.matshow(d)
    w = weight(distance)
    # fig, ax = plt.subplots(ncols=2)
    # ax[1].matshow(feature)
    # X, Y = window[1] .+ min[1], window[2] .+ min[2]
    # ax[2].matshow(target[X, Y])
    # Convert from surface coordinates to feature indices
    offset = first.(window) .+ (length.(window) .÷ 2)
    return min[1] + offset[1], min[2] + offset[2], w
end

function featureoverlap2d_ensemble(feature, target, searchwidth, N::Number, members=1, metric=χ²)
    minima = []
    perturbed = similar(feature)
    #for member in ProgressBar(1:members)
    window, _ = shiftsinrange(feature, target, searchwidth)
    metric = L2
    for member in 1:members
        # Perturb the target to add some randomness
        perturbate!(perturbed, feature)
        perturbed = imfilter(perturbed, Kernel.gaussian(3))
        perturbed ./= N
        distance = slidingwindow(perturbed, target, feature.indices, metric, searchwidth)
        min = argmin(distance)
        fig, ax = plt.subplots()
        ax.matshow(distance)
        ax.scatter(min[2], min[1])
        fig, ax = plt.subplots(ncols=2)
        ax[1].matshow(perturbed)
        X, Y = window[1] .+ min[1], window[2] .+ min[2]
        ax[2].matshow(target[X, Y])
        push!(minima, min)
    end
    println(minima[1])
    # Convert from surface coordinates to feature indices
    return minima
    # Transform list of ranges to matrix of midpoints
    # mid = [m .|> middle for m in minima]
    # mid = transpose(hcat(mid...))
    # fig, ax = subplots()
    # ax.scatter(mid[:, 1], mid[:, 2], color="r", alpha=0.5)
    # midpoint = round.(Int, median(mid, dims=1))
    # ax.scatter(midpoint[1], midpoint[2], color="g")
    # x, y = feature.indices
    # wx = (x[end] - x[1])÷2
    # wy = (y[end] - y[1])÷2
    # mid = midpoint
    # [(mid[1]-wx):(mid[1]+wx), (mid[2]-wy):(mid[2]+wy)]
end

function squarefrommidpoint(midpoint, width)
    squarefrommidpoint(midpoint, width...)
end

function squarefrommidpoint(midpoint, xwidth, ywidth)
    wx = xwidth ÷ 2
    wy = ywidth ÷ 2
    [(midpoint[1]-wx):(midpoint[1]+wx), (midpoint[2]-wy):(midpoint[2]+wy)]
end

"""
    slidingwindow()

Compares `reference` to sliding windows of equal size of `target`
across a square of width `searchwidth` using the metric `metric`.
The search window is centered at the center of reference, under the
assumption that the reference is a view of a larger matrix.
Returns the resulting potential surface.
"""
function slidingwindow(reference::SubArray, target::AbstractMatrix, metric::Function, searchwidth::Int)
    return slidingwindow(reference, target, reference.indices, metric, searchwidth)
end

function slidingwindow(reference::AbstractMatrix, target::AbstractMatrix, 
                       indices, metric::Function, searchwidth::Int)
    window, shifts = shiftsinrange(reference, target, indices, searchwidth)
    distance = Matrix(undef, shifts[1], shifts[2])
    for i in 1:shifts[1], j in 1:shifts[2]
        # The .- 1 accounts for zero shift
        distance[i, j] = metric(reference, target[window[1] .+ i .- 1,
                                                  window[2] .+ j .- 1])
    end
    return distance
end


function shiftsinrange(template::SubArray, searchregion::AbstractMatrix, searchwidth::Int)
    shiftsinrange(template, searchregion, template.indices, searchwidth)
end

function shiftsinrange(template, searchregion::AbstractMatrix, indices, searchwidth::Int)
    @assert searchwidth < size(template, 1) "$searchwidth > $(size(template, 1))"
    @assert searchwidth < size(template, 2) "$searchwidth > $(size(template, 2))"
    # Handle boundary issues by clipping the search range to within [1, size target]
    start = clip.(first.(indices) .- searchwidth, size(searchregion))
    stop  = clip.(last.(indices)  .+ searchwidth, size(searchregion))
    # Construct the indices defining the initial window in the upper left
    # corner of the search box
    window = range.(start, start .+ size(template) .- 1, step=1)
    # The number of shifts of the window possible within the search range
    shifts = stop .- start .- size(template)
    return window, shifts
end

# function featureoverlap_perturbed(feature, target, searchwidth, indices)
    # x, y = indices 
    # widthx, widthy = size(feature)
    # startx  = clip(first(x) - searchwidth, size(target, 1))
    # stopx   = clip(last(x)  + searchwidth, size(target, 1))
    # starty  = clip(first(y) - searchwidth, size(target, 2))
    # stopy   = clip(last(y)  + searchwidth, size(target, 2))
    # windowx = startx:(startx + widthx - 1)
    # windowy = starty:(starty + widthy - 1)
    # lengthx, lengthy = (stopx - startx - widthx,
                        # stopy - starty - widthy)
    # overlap = zeros((lengthx, lengthy))
    # for j in 1:lengthy, i in 1:lengthx
        # shiftx = i - 1  # Account for zero shift
        # shifty = j - 1
        # diff = (feature .- target[windowx .+ shiftx,
                                  # windowy .+ shifty]).^2 ./ feature
        # overlap[i, j] = sum(diff)
    # end
    # #return overlap

    # # Find the best overlap
    # indices = argmin(overlap)
    # x = range(startx + indices[1] + 1, length=widthx)
    # y = range(starty + indices[2] + 1, length=widthy)
    # [x, y]
# end


function perturbate(raw::AbstractArray)
    perturbed = similar(raw)
    for i in eachindex(raw)
        mu = raw[i]# < 1 ? 1 : mat[i]
        perturbed[i] = rand(Poisson(mu))
    end
    perturbed
end

function perturbate!(perturbed::AbstractArray, raw::AbstractArray)
    for i in eachindex(raw)
        mu = raw[i]# < 1 ? 1 : mat[i]
        #perturbed[i] = mu < 2 ? 0 : rand(Poisson(mu))
        perturbed[i] = rand(Poisson(mu))
    end
end


function arraybox(array, startx::T, stopx::T, starty::T, stopy::T) where T <: Integer
    subarray = view(array, startx:stopx, starty:stopy)
end

# function splitregions2d(array::OMatrix{T}, numberofregions; width=nothing, limit=50) where T
#     splitregions2d(array, numberofregions, width=width, limit=limit)
# end

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

function downsample(mat::AbstractMatrix{T}, blocksize=2) where T
    xlen, ylen = size(mat)
    newx, newy = xlen÷blocksize, ylen÷blocksize
    @assert(xlen%blocksize == 0, "Dimension 1 is not divisible by blocksize, $xlen")
    @assert(ylen%blocksize == 0, "Dimension 2 is not divisible by blocksize, $ylen")
    new = Matrix{Float64}(undef, newx, newy)
    s = 0
    for i in 1:blocksize:xlen, j in 1:blocksize:ylen
        s = 0
        for i_ in i:(i+blocksize-1), j_ in j:(j+blocksize-1)
            s += mat[i_, j_]
        end
        new[i÷blocksize + 1, j÷blocksize + 1] = 1/blocksize*s
    end
    return new
end

function neighbourreduce(mat::AbstractMatrix{T}; width=3, threshold=10) where T
    @assert isodd(width)
    xlen, ylen = size(mat)
    new = Matrix{Float64}(undef, ylen, xlen)
    s = 0
    w = width÷2
    for i in 1:xlen, j in 1:ylen
        s = 0
        for i_ in (i-w):(i+w), j_ in (j-w):(j+w)
            if 0 < i_ < xlen && 0 < j_ < ylen
                s += mat[i_, j_]
            end
        end
        if s ≥ threshold
            new[i, j] = mat[i, j]
        else
            new[i, j] = 0.0
        end
    end
    return new
end

function annotate_plot_2d(mat, x, y, w=nothing; ax=nothing)
    if ax ≡ nothing
        fig, ax = plt.subplots()
    end

    ax.matshow(transpose(mat) .|> log10)
    ax.scatter(x, y, marker="x", c="r")
    if w !== nothing
        for (i, w_) in enumerate(w)
            ax.annotate(string(w_), (x[i], y[i]))
        end
    end
    ax.invert_yaxis()
end

function compare_2d(x, y;kernel=Kernel.gaussian(3), metric::Function=L2)
    x = imfilter(x, kernel)
    y = imfilter(y, kernel)
    # Normalize y to x
    # y .= y.*sum(x)/sum(y)
    # Normalize each
    x ./= sum(x)
    y ./= sum(y)
    return metric(x, y)
    #return (x .- y).^2
end

function compare_methods(reference, target; neighbourwidth=9, neighbourthreshold=8)
    edges = histogram_edges(reference[!, :e], reference[!, :Δe], nbins=1000)
    refmat_raw = JSort.histogram(reference[!, :e], reference[!, :Δe], edges...)
    tarmat_raw = JSort.histogram(target[!, :e], target[!, :Δe], edges...)

    Reds = alphamap("Reds")
    Blues = alphamap("Blues")

    fig, ax = plt.subplots()
    ax.set_title("Original")
    ax.pcolormesh(refmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                  cmap=Reds, alpha=0.7)
    ax.pcolormesh(tarmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                  cmap=Blues, alpha=0.7)
    ax.set_ylim([0, 600])
    ax.set_xlim([0, 800])

    refmat = neighbourreduce(refmat_raw, width=neighbourwidth, threshold=neighbourthreshold)
    tarmat = neighbourreduce(tarmat_raw, width=neighbourwidth, threshold=neighbourthreshold)
    points = featurealign2d_weighted(refmat, tarmat, numregions=15, width=nothing,
                                 limit=1500, searchwidth=40, weight=ThresholdSum(.2))
    rpx, tpx, rpy, tpy, weights = points
    annotate_plot_2d(tarmat, tpx, tpy, weights)
    edgesx, edgesy = edges
    rpx, tpx = edgesx[rpx], edgesx[tpx]
    rpy, tpy = edgesy[rpy], edgesy[tpy]

    lm = pyimport("sklearn.linear_model")

    models = [leastsquares, 
              (x, y) -> begin
                lr = lm.TheilSenRegressor(fit_intercept=true)
                X = reshape(x, :, 1)
                lr.fit(X, y)
                return [lr.intercept_; lr.coef_]
              end,
              (x, y) -> begin
                  lr = lm.RANSACRegressor()
                  X = reshape(x, :, 1)
                  lr.fit(X, y)
                  return [lr.estimator_.intercept_; lr.estimator_.coef_]
              end]
    names = ["OLS", "TheilSen", "RANSAC"]
    metric = (x, y) -> splitmetric(x, y, metric=L2)

    for (model, name) in zip(models, names)
        coeffx = model(tpx, rpx)
        coeffy = model(tpy, rpy)

        px = Poly(coeffx)
        py = Poly(coeffy)

        e  = px.(target[!, :e])
        de = py.(target[!, :Δe])

        calibmat = JSort.histogram(e, de, edges...);
        err1 = compare_2d(refmat, calibmat, metric=metric)
        err2 = compare_2d(refmat, tarmat, metric=metric)
        println("$name - Error ratio: $(err1/err2)")

        fig, ax = plt.subplots(ncols=2, nrows=2)
        plotresiduals(tpx, rpx, coeffx, ax=ax[1:2])
        plotresiduals(tpy, rpy, coeffy, ax=ax[3:4])
        fig.suptitle(name)

        fig, ax = plt.subplots(figsize=(15, 15))
        ax.set_title("Aligned")
        ax.pcolormesh(refmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Reds, alpha=0.7)
        ax.pcolormesh(calibmat', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Blues, alpha=0.7)
        axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
        axins.pcolormesh(refmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Reds, alpha=0.7)
        axins.pcolormesh(calibmat', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Blues, alpha=0.7)
        axins.set_xlim(550, 800)
        axins.set_ylim(80, 220)
        ax.indicate_inset_zoom(axins)
        axins = ax.inset_axes([0.0, 0.5, 0.47, 0.47])
        axins.pcolormesh(refmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Reds, alpha=0.7)
        axins.pcolormesh(calibmat', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Blues, alpha=0.7)
        axins.set_xlim(0, 200)
        axins.set_ylim(400, 600)
        ax.indicate_inset_zoom(axins)
    end

    loss = ["absolute_loss", "squared_loss"]
    models = [(x, y, w) -> begin
                  X = reshape(x, :, 1)
                  lr.fit(X, y, w)
                  return [lr.estimator_.intercept_; lr.estimator_.coef_]
              end
                for lr in 
              [lm.RANSACRegressor(loss=l) for l in loss]]
    names = ["RANSAC "*l for l in loss]

    # Weighted
    for (name, model) in zip(names, models)
        coeffx = model(tpx, rpx, weights)
        coeffy = model(tpy, rpy, weights)

        px = Poly(coeffx)
        py = Poly(coeffy)

        e  = px.(target[!, :e])
        de = py.(target[!, :Δe])

        calibmat = JSort.histogram(e, de, edges...);
        err1 = compare_2d(refmat, calibmat, metric=metric)
        err2 = compare_2d(refmat, tarmat, metric=metric)
        println("$name - Error ratio: $(err1/err2)")

        fig, ax = plt.subplots(ncols=2, nrows=2)
        plotresiduals(tpx, rpx, coeffx, ax=ax[1:2])
        plotresiduals(tpy, rpy, coeffy, ax=ax[3:4])
        fig.suptitle(name)

        fig, ax = plt.subplots(figsize=(15, 15))
        ax.set_title("Aligned")
        ax.pcolormesh(refmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Reds, alpha=0.7)
        ax.pcolormesh(calibmat', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Blues, alpha=0.7)
        axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])
        axins.pcolormesh(refmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Reds, alpha=0.7)
        axins.pcolormesh(calibmat', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Blues, alpha=0.7)
        axins.set_xlim(550, 800)
        axins.set_ylim(80, 220)
        ax.indicate_inset_zoom(axins)
        axins = ax.inset_axes([0.0, 0.5, 0.47, 0.47])
        axins.pcolormesh(refmat_raw', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Reds, alpha=0.7)
        axins.pcolormesh(calibmat', norm=PyPlot.matplotlib.colors.LogNorm(),
                      cmap=Blues, alpha=0.7)
        axins.set_xlim(0, 200)
        axins.set_ylim(400, 600)
        ax.indicate_inset_zoom(axins)
    end
    println("bump")
end

function alphamap(name::String)
PyCall.py"""
import numpy as np
import matplotlib.pylab as pl
from matplotlib.colors import ListedColormap

def alpha_map(name):
    # Choose colormap
    cmap = pl.cm.get_cmap(name)

    # Get the colormap colors
    my_cmap = cmap(np.arange(cmap.N))

    # Set alpha
    my_cmap[:,-1] = np.linspace(0, 1, cmap.N)

    # Create new colormap
    my_cmap = ListedColormap(my_cmap)
    return my_cmap
"""
return PyCall.py"alpha_map"(name)
end
