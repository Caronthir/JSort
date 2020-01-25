using JSort
using Statistics
using PyCall
# @pyimport scipy.signal as signal


function findpeaks(spectrum; numpeaks=2)::Vector{Int64}
    signal = pyimport("scipy.signal")
    spectrum_norm = spectrum/mean(spectrum)
    prominence = 20
    peaks = []
    while length(peaks) < numpeaks
        peaks, prop = signal.find_peaks(spectrum_norm, prominence=prominence)
        prominence -= 0.1
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
                      lowsearchwidth=50, highsearchwidth=500, referenceindex=1, order=2)
    if isnothing(lowregion)
        lowregion = 1:length(spectra[1])
    end

    reference = spectra[referenceindex]
    coefficients = Vector{Float64}[]
    aligned = Vector{Int}[]
    for i in eachindex(spectra)
        if i == referenceindex
            push!(coefficients, [0.0, 1.0, 0.0])
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
        push!(aligned, gainshift(target, shiftgain...))

        # Convert from index basis to energy basis
        shift = first(energyedges)
        gain = step(energyedges)
        targetpeaks = targetpeaks.*gain .+ shift
        referencepeaks = referencepeaks.*gain .+ shift
        shiftgain = leastsquares(targetpeaks, referencepeaks, order=order)
        push!(coefficients, shiftgain)
    end
    coefficients, aligned
end


function splitregion(region, numberofregions; width=150)::Vector{Tuple{Int64, Int64}}
    # Create [start, stop] indices for each region
    step = floor.(Int, length(region)/numberofregions)
    indices = []
    for i in 1:numberofregions
        start = (i-1)*step + 1
        stop = i*step + 1
        middle = start + (stop-start)/2
        lower = middle-width/2 |> x -> floor(Int64, x)
        upper = middle+width/2 |> x -> floor(Int64, x)

        (lower < 1) && (lower = 1)
        (upper > length(region)) && (upper = length(region))

        push!(indices, (lower, upper))
    end
    #indices = [((i-1)*step + 1, i*step + 1) for i in 1:numberofregions-1]
    #push!(indices, ((numberofregions-1)*step+1, length(region)))
    return indices
end


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
        lower = clip(start - searchwidth, length(target))
        upper = clip(stop + searchwidth, length(target))
        feature_target = featureshift(feature, target, start=lower,
                                      stop=upper)
        push!(referencepeaks, middle(start:stop))
        push!(targetpeaks, middle(feature_target))
    end
    (referencepeaks, targetpeaks) .|> x -> x .+ first(roi) .- 1
end

function clip(x::T, upper::T)::T where T
    if x < 1
        one(x)
    elseif x > upper
        upper
    else
        x
    end
end

middle(x) = length(x)/2 + first(x) |> x -> round(Int64, x)
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
    range(start + index + 1, length=width)
end
