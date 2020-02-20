using .JSort
using PyCall
import PyPlot; const plt = PyPlot;


function localalign(detector::ChannelDetector,
                    region; nbins=5000, emin=-150, emax=15000,
                    kwargs...)
    full, edges = histogram(detector, nbins=nbins)
    chns = goodchannels(detector)

    signal = pyimport("scipy.signal")
    smoother(x) = signal.savgol_filter(x, 51, 4)

    eclip(x) = clip(x, emin, emax)

    kw = get(kwargs, :a1, Dict())
    if get(kw, :plot, false)
        plothist(full, edges)
    end
    delta = get(kw, :delta, [3500, 3500])
    searchregion = eclip(region[1]-delta[1]):eclip(region[2]+delta[2])

    println("Aligning")

    coeff, _ = alignspectra(full, edges, lowregion=searchregion, lowsmoother=smoother,
                            lownumregions=get(kw, :nregions, 2),
                            lowregionwidth=get(kw, :swidth, 200),
                            lowsearchwidth=get(kw, :lwidth, 60),
                            plot=get(kw, :plot, false), order=1)
    coefficients = repeat([0.0 1.0], outer=(33, 1))
    coefficients[chns, :] .= hcat(coeff...)'
    calib = GammaCalibrator(coefficients)
    aligned = calib(detector, inplace=false)

    if get(kw, :plot, false)
        plot(aligned)
    end

    alignedhist = histogram(aligned, edges)

    kw = get(kwargs, :a2, Dict())
    delta = get(kw, :delta, [3500, 3500])
    searchregion = eclip(region[1]-delta[1]):eclip(region[2]+delta[2])

    coeff, _ = alignspectra(alignedhist, edges, lowregion=searchregion, lowsmoother=smoother,
                            lownumregions=get(kw, :nregions, 3),
                            lowregionwidth=get(kw, :swidth, 200),
                            lowsearchwidth=get(kw, :lwidth, -10),
                            plot=get(kw, :plot, false), order=1)

    coefficients = repeat([0.0 1.0], outer=(33, 1))
    coefficients[chns, :] .= hcat(coeff...)'
    combine!(calib, GammaCalibrator(coefficients))
    aligned = calib(detector, inplace=false)

    if get(kw, :plot, false)
      plot(aligned)
    end

    #S = sum(aligned)
    #hist = histogram(S, edges)
    #plot(S)
    # Energy - index is wrong
    #peak = findpeaks(hist[region], numpeaks=1)[1]
    #peak += first(region)
    #@show peak
    #@show peaks = map(x -> x(peak), invcalib.poly)
    #ax.scatter(peaks[chns], 1.5:length(chns)+0.5, marker="x", color="r", s=2)

    return calib
end


function regionalign(d::ChannelDetector, regions; nbins=1000,
                     emin=100, emax=2500, kwargs...)
    strips, edges = JSort.histogram(d, nbins=nbins)
    shift, gain = first(edges), step(edges)
    indextoe(x) = shift + gain*x
    etoindex(x) = round(Int, -shift/gain + 1/gain*x)
    #regions = [map(etoindex, region) for region in regions]
    eclip(x) = clip(x, emin, emax)
    delta = 400

    peaks = []
    for region in regions
        low, high = region[1], region[2]
        # Align all strips to the reference strip
        lowregion = eclip(low-delta):eclip(high+delta)
        @show lowregion
        println("Aligning")
        al, calib, _ = align(d, lowregion=lowregion, lownumregions=2,
                             lowregionwidth=200,
                             lowsearchwidth=80, order=1, kwargs...);
        println("Done aliging")
        println("Region peaking")
        # Find the peak(s) in the aligned region
        peak = regionpeaks(al, [region])
        println("Done region peaking")
        # Backshift the peaks to the original coordinate system
        peak = backshift(peak, goodchannels(d), invert(calib))
        push!(peaks, peak)
    end
    peaks
end

function regionpeaks(d::ChannelDetector, regions; nbins=5000)
    strips, edges = JSort.histogram(d, nbins=nbins)
    shift, gain = first(edges), step(edges)
    indextoe(x) = shift + gain*x
    etoindex(x) = round(Int, -shift/gain + 1/gain*x)
    regions = [map(etoindex, region) for region in regions]
    println("Regioning for $regions")
    peaks = regionpeaks(strips, regions)
    [map(indextoe, peak) for peak in peaks]
end

function regionpeaks(strips, regions)
    peaks = Vector{Int}[]
    for strip in strips
        peaks_ = regionpeaks(strip, regions)
        push!(peaks, peaks_)
    end
    peaks
end

function regionpeaks(strip::AbstractVector{<:Number}, regions)::Vector{Int}
    peaks = Int[]
    for region in regions
        peak = regionpeak(strip, region[1], region[2])
        push!(peaks, peak)
    end
    return peaks
end

function regionpeak(strip::AbstractVector{<:Number}, lower::Int, upper::Int)::Int
    region = strip[lower:upper]
    peak = findpeaks(region, numpeaks=1)[1]
    return peak + lower - 2 #BUG -2 is indicative of a offset bug
end

function backshift(points, channels, calib)
    shiftedpoints = similar(points)
    for i in 1:length(points)
        shiftedpoints[i] = calib[channels[i]](points[i])
    end
    shiftedpoints
end

function backshift(points::AbstractVector, coefficients::AbstractVector)
    shiftedpoints = similar(points)
    for i in 1:length(points)
        shift, gain = coefficients[i]
        poly(x) = x/gain - shift/gain
        shiftedpoints[i] = poly(points[i])
    end
    shiftedpoints
end

function backshift(point::Number, coefficients::AbstractVector)
    shiftedpoints = []
    for i in 1:length(points)
        shift, gain = coefficients[i]
        poly(x) = x/gain - shift/gain
        push!(shiftedpoints, poly(point))
    end
    shiftedpoints
end

function clip(x, lower, upper)
    if x < lower
        lower
    elseif x > upper
        upper
    else
        x
    end
end

function plothist(strips, edges; points=nothing)
    mat = cat(strips'..., dims=1)
    fig, ax = plt.subplots()
    ax.pcolormesh(edges, 1:length(strips)+1, mat .|> log10)
    if !isnothing(points)
        for points_ in points
            for (i, point) in enumerate(points_)
                ax.scatter(point, i+0.5, marker="x", color="r", s=3)
            end
        end
    end
    ax
end

function plothist(det::ChannelDetector; nbins=1000, points=nothing)
    det_strips, det_edges = JSort.histogram(det, nbins=nbins)
    plothist(det_strips, det_edges, points=points)
end
