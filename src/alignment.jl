using .JSort
using Statistics
using PyCall
using DocStringExtensions
# @pyimport scipy.signal as signal


function findpeaks(spectrum; numpeaks=2)::Vector{Int}
    signal = pyimport("scipy.signal")
    spectrum_norm = spectrum/mean(spectrum)
    prominence = 20
    peaks = []
    i = 0
    prev = 0
    Δ = 0.1
    while length(peaks) != numpeaks
        peaks, prop = signal.find_peaks(spectrum_norm, prominence=prominence)
        if length(peaks) < numpeaks
            prominence -= Δ
            if prev > numpeaks && i > 5
                Δ /= 10
            end
        else
            prominence += Δ
            if prev < numpeaks && i > 5
                Δ /= 10
            end
        end

        i += 1
        if i > 2000
            break
        end
        prev = length(peaks)
    end
    return peaks[1:numpeaks] .+ 1
end

function findpeaks(spectrum, region; numpeaks=2)::Vector{Int64}
    # TODO: Can be made more efficient with views
    # TODO: Maybe support discrete regions
    signal = pyimport("scipy.signal")
    spectrum_norm = spectrum[region]/mean(spectrum[region])
    prominence = 20
    peaks = []
    while length(peaks) < numpeaks
        peaks, prop = signal.find_peaks(spectrum_norm, prominence=prominence)
        prominence -= 0.1
    end
    x = Int[]
    for peak in peaks[1:numpeaks]
        push!(x, peak + first(region))
    end
    return x .+ 1  # For indexing difference Python/Julia
end

function lastnonzero(spectrum)
    for (i, elem) in enumerate(reverse(spectrum))
        if elem ≠ 0
            return i
        end
    end
end

function goodchannels(matrix::T; lastzerolimit=100, deviationlimit=10)::Vector{Int} where T <: Dict
    good = Int[]
    means = Float64[]
    sizes = Int[]
    for (channel, values) in matrix
        last = lastnonzero(values)
        isnothing(last) && continue
        last > lastzerolimit && continue
        push!(good, channel)
        push!(means, mean(values))
        push!(sizes, size(values, 1))
    end
    # Second pass to remove those with abnormal mean
    M = mean(means)
    deviation = @. 1/(sizes - 1)*√((means - M)^2)

    good[deviation .< deviationlimit] |> sort
end

function goodchannels(matrix; lastzerolimit=100, deviationlimit=10)::Vector{Int}
    good = Int[]
    means = Float64[]
    for i in 1:size(matrix, 2)
        last = lastnonzero(matrix[:, i])
        isnothing(last) && continue
        last > lastzerolimit && continue
        push!(good, i)
        push!(means, mean(matrix[:, i]))
    end
    # Second pass to remove those with abnormal mean
    M = mean(means)
    deviation = 1/(size(matrix, 2) -1)*.√((means .- mean(means)).^2)
    good[deviation .< deviationlimit]
end


function linearalign(X, ref; numpeaks=2)
    peaks = findpeaks(X, numpeaks=numpeaks)
    refpeaks = findpeaks(ref, numpeaks=numpeaks)
    leastsquares(peaks, refpeaks)
end

function linearalign(X, ref, region; numpeaks=2)
    peaks = findpeaks(X, region, numpeaks=numpeaks)
    refpeaks = findpeaks(ref, region, numpeaks=numpeaks)
    leastsquares(peaks, refpeaks)
end


function alignspectra(spectra, energyedges; lowregion=nothing, highregion=nothing,
                      lowregionwidth=250, highregionwidth=120, lownumregions=10,
                      highnumregions=2, lowsmoother=nothing, highsmoother=nothing,
                      lowsearchwidth=200, highsearchwidth=500, referenceindex=1, order=2,
                      plot=false)
    shift = first(energyedges)
    gain = step(energyedges)
    indextoe = Poly([shift, gain])
    etoindex = Poly([-shift/gain, 1/gain])
    etoindexint = Int∘round∘etoindex

    @show lowsearchwidth
    # Convert from energy basis to index basis
    lowregionwidth = etoindexint(lowregionwidth)
    highregionwidth = etoindexint(highregionwidth)
    lowsearchwidth = etoindexint(lowsearchwidth)
    highsearchwidth = etoindexint(highsearchwidth)

    @show lowsearchwidth
    if lowsearchwidth < 180
        lowsearchwidth = 50
    end

    if isnothing(lowregion)
        lowregion = 1:length(spectra[1])
    else
        lowregion = etoindexint(first(lowregion)):etoindexint(last(lowregion))
        if plot
          plt.axvline(x=first(lowregion) |> indextoe, c="w")
          plt.axvline(x=last(lowregion) |> indextoe, c="w")
        end
    end

    if !isnothing(highregion)
        highregion = etoindexint(first(highregion)):etoindexint(last(highregion))
        if plot
          plt.axvline(x=minimum(highregion) |> indextoe, c="w", linestyle="--")
          plt.axvline(x=maximum(highregion) |> indextoe, c="w", linestyle="--")
        end
    end

    reference = spectra[referenceindex]
    coefficients = Vector{Float64}[]
    coefficients⁻¹ = Vector{Float64}[]
    aligned = Vector{Int}[]
    default = if order == 0
        [0.0]
    elseif order == 1
        [0.0, 1.0]
    else
        [0.0, 1.0, 0.0]
    end
    flag = false
    for i in eachindex(spectra)
        if i == referenceindex
            push!(coefficients, default)
            push!(coefficients⁻¹, default)
            push!(aligned, reference)
            continue
        end

        target = spectra[i]
        referencepeaks, targetpeaks = featurealign(reference, target,
                                                   width=lowregionwidth,
                                                   searchwidth=lowsearchwidth,
                                                   roi=lowregion, numregions=lownumregions,
                                                   smoother=lowsmoother)

        if !isnothing(highregion)
            refhigh, tarhigh = featurealign(reference, target,
                                            width=highregionwidth,
                                            searchwidth=highsearchwidth,
                                            roi=highregion, numregions=highnumregions,
                                            smoother=highsmoother)
            push!(referencepeaks, refhigh...)
            push!(targetpeaks, tarhigh...)
        end
        shiftgain = leastsquares(targetpeaks, referencepeaks, order=order)
        #push!(aligned, gainshift(target, shiftgain...))

        # Convert from index basis to energy basis
        # targetpeaks    = targetpeaks.*gain .+ shift
        # referencepeaks = referencepeaks.*gain .+ shift
        targetpeaks    = indextoe.(targetpeaks)
        referencepeaks = indextoe.(referencepeaks)
        if plot
          plt.scatter(targetpeaks, repeat([i+0.5], length(targetpeaks)), color="r", marker="v")
          if !flag
              plt.scatter(referencepeaks, repeat([referenceindex+0.5], length(referencepeaks)), color="r", marker="v")
              flag = true
          end
        end
        shiftgain = leastsquares(targetpeaks, referencepeaks, order=order)
        shiftgain⁻¹ = leastsquares(referencepeaks, targetpeaks, order=order)
        push!(coefficients, shiftgain)
        push!(coefficients⁻¹, shiftgain⁻¹)
    end
    coefficients, coefficients⁻¹
end


"""

$(SIGNATURES)

The array `region` is split into `numberofregions` regions, each with a length of
`width` measured from the middle of each sub-region. The subregions can overlap, but are
clipped if they fall outside the range of the array.

    |...............|
           ↓ Split into `numberofregions` = 3
    |.....|.....|.....|
           ↓ From the middle of each subregion, select a subsubregion of `width` = 3
    |.¦...¦. | .¦...¦.| .¦...¦.|
           ↓ Return the (start, stop) of each subsubregion
    [(2, 4), (7, 9), (12, 14)]
# Examples
```julia-repl
julia> splitregion(1:15, 3, width=3)
3-element Array{Tuple{Int64,Int64},1}:
 (2, 4)
 (7, 9)
 (12, 14)
```
"""
function splitregion(region, numberofregions::Int; width=150)::Vector{Tuple{Int64, Int64}}
    # Create [start, stop] indices for each region
    step = floor.(Int, length(region)/numberofregions)
    indices = []
    for i in 1:numberofregions
        start = (i-1)*step + 1
        stop = i*step + 1
        middle = start + (stop-start)/2
        lower = middle-width/2 |> x -> floor(Int64, x)
        upper = middle+width/2 |> x -> floor(Int64, x) - 1

        (lower < 1) && (lower = 1)
        (upper > length(region)) && (upper = length(region))

        push!(indices, (lower, upper))
    end
    #indices = [((i-1)*step + 1, i*step + 1) for i in 1:numberofregions-1]
    #push!(indices, ((numberofregions-1)*step+1, length(region)))
    return indices
end


"""

$(SIGNATURES)

The arrays `reference` and `target` are compared in the region of interest `roi`.
This region is split into `numregions` subregions, where each subregion is `width` long.
A subregion is cut out of the `reference` array and compared to the `target` along a window
which is +- `searchwidth` longer than the width of the `feature`. The corresponding region
is found in the `target`. The middle point of point the reference feature and the target
feature is returned for each region.

"""
function featurealign(reference, target; roi=nothing, numregions=10, width=50,
                      searchwidth=50, smoother=nothing)
    if isnothing(roi)
        roi = 1:length(reference)
    end
    if !isnothing(smoother)
        reference = smoother(reference)
        target = smoother(target)
    end
    reference = reference[roi]
    target = target[roi]

    regions = splitregion(reference, numregions, width=width)
    referencepeaks = Int64[]
    targetpeaks = Int64[]
    for (start, stop) in regions
        feature = reference[start:stop]
        # Find the start and stop indices in the target whence to compare the feature
        lower = clip(start - searchwidth, length(target))
        upper = clip(stop  + searchwidth, length(target))
        feature_target = featureshift(feature, target, start=lower,
                                      stop=upper)
        push!(referencepeaks, middle(start:stop))
        push!(targetpeaks, middle(feature_target))
    end
    (referencepeaks, targetpeaks) .|> x -> x .+ first(roi) .- 1
end

function clip(x, upper)
    if x < 1
        one(x)
    elseif x > upper
        upper
    else
        x
    end
end

middle(x) = length(x)/2 + first(x) |> x -> round(Int64, x)
"""
$(SIGNATURES)

Compares the `feature` along the `target` from the indices `start` to `stop`.
A window of equal size as the `feature` is created and slid along the entire
start:stop range with the error computed as the MSE. The point of minimum MSE
is selected.
"""
function featureshift(feature, target; start=nothing, stop=nothing)
    isnothing(start) && (start = 1)
    isnothing(stop) && (stop = length(target))

    width = length(feature)
    window = start:(start + width - 1)
    errors = zeros(stop - start - width)
    for i in eachindex(errors)
        shift = i - 1  # For the zero shift
        diff = (feature - target[window.+shift]).^2
        errors[i] = sum(diff./feature.^2)
    end
    # TODO Use peakdetection?
    index = argmin(errors)
    #index = findpeaks(1/errors, numpeaks=1)[1]
    range(start + index, length=width)
end
