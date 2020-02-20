using .JSort
using PyCall
import PyPlot; const plt = PyPlot;

function iterativealign(parameters::Parameters, cascades::AbstractVector{GammaCascade};
                        kwargs...)
    ungated = QuickReader(parameters);
    gated = QuickReader(parameters);
    iterativealign(ungated, gated, cascades; kwargs...)
end

function iterativealign(ungated::QuickReader, gated::QuickReader,
                        cascades::AbstractVector{GammaCascade};
                        nbins=5000, nbins2=3000, kwargs...)
    # TODO This is complicated with a lot of parts. Think more carefully
    # For each γ cascade, gate on the corresponding excitation energy window
    # And fit each γ decay
    # TODO Needs to perform some initial calibration so that the channel
    #      numbers aren't random
    default = GammaCalibrator([-913.5, 5.6])
    println("Reading ungated")
    read(ungated)
    ungateddetector = ChannelDetector(ungated)
    default(ungateddetector)
    signal = pyimport("scipy.signal")
    smoother(x) = signal.savgol_filter(x, 51, 4)
    chns = goodchannels(ungateddetector)
    # ungatedhist, edges = histogram(ungateddetector, nbins=nbins)
    # Contains the peaks for all gamma cascades for all spectra
    allpeaks = Vector{Vector{Float64}}[]
    for cascade in cascades
        println("Reading gated for $(cascade.gate)")
        read(gated, cascade.gate)
        gateddetector = ChannelDetector(gated)
        default(gateddetector)
        # gatedhist = histogram(gateddetector, edges)
        # Contains the peaks for all decays in a cascade
        # for all spectra
        peaks = Vector{Float64}[]
        individual = true
        for (e, region) in zip(cascade.e, cascade.regions)
            nbins = get(cascade.keywords, :bins1, 5000)
            nbins2 = get(cascade.keywords, :bins2, 3000)

            region = first(region):last(region)
            calib = localalign(ungateddetector, region; nbins=nbins, cascade.keywords[:align]...)
            aligned = calib(gateddetector, inplace=false)
            ax = plot(aligned, nbins=nbins)
            @show region
            plt.axvline(x=first(region), c="r")
            plt.axvline(x=last(region), c="r")

            total = sum(aligned)
            plot(total, nbins=nbins2)
            plt.axvline(x=first(region), c="k")
            plt.axvline(x=last(region), c="k")

            inv = invert(calib)
            peaks_ = if !individual
                peaks_ = localpeak(total, region, nbins=nbins2, smoother=smoother)
                peaks_ = inv(peaks_)[chns]
            else
                peaks_ = []
                for (detector, chn) in zip(aligned.data, 1:33)
                    chn ∉ chns && continue
                    peak = localpeak(detector, region, nbins=nbins2, smoother=smoother)
                    peak = inv[chn](peak)
                    push!(peaks_, peak)
                end
                peaks_
            end

            ax = plot(ungateddetector)
            ax.scatter(peaks_, 1.5:length(peaks_)+1, marker="x", color="r")

            push!(peaks, peaks_)
        end
        push!(allpeaks, peaks)
    end
    ax = plot(ungateddetector)
    for peaks in allpeaks
        for peaks_ in peaks
            ax.scatter(peaks_, 1.5:length(peaks_)+1, marker="x", color="r")
        end
    end
    allpeaks
end

function localpeak(detector::GammaDetector, region; nbins=3000, smoother=nothing)
    hist, edges = histogram(detector, nbins=nbins)
    index2e, e2index, e2indexint = transformations(edges)

    if !isnothing(smoother)
        hist = smoother(hist)
    end

    regioni = e2indexint(first(region)):e2indexint(last(region))
    peaki = findpeaks(hist[regioni], numpeaks=1)[1]
    #@show hist[regioni][peaki]
    peaki += first(regioni) |> Int
    #@show hist[peaki]
    #fig, ax = plt.subplots()
    #ax.step(edges[regioni], hist[regioni])
    #ax.plot(edges[regioni], smoothed[regioni])
    #ax.scatter(edges[peaki], hist[peaki])
    #plt.step(edges[regioni], hist[regioni])
    #plt.scatter(edges[peaki], hist[peaki])
    peak = index2e(peaki)
    #plt.scatter(peak, hist[peaki], marker="v", color="r", s=50)
    peak
end

function iterativealign(cascade::GammaCascade, ungatedhist,
                        edges, gated::QuickReader, peakwidth=50
                        ;kwargs...)
    println("Reading gated for $(cascade.gate)")
    read(gated, gate=cascade.gate)
    gateddetector = ChannelDetector(gated)
    gatedhist = histogram(gateddetector, edges)
    # Contains the peaks for all decays in a cascade
    # for all spectra
    peaks = Float64[]
    for e in cascade.e
        region = (e-peakwidth, e+peakwidth)
        peak = regionalign(ungatedhist, gatedhist,
                           edges, region;
                           kwargs...)
        push!(peaks, peak)
    end
    peaks
end

function transformations(range::AbstractRange)
    shift = first(range)
    gain = step(range)
    indextoe = Poly([shift, gain])
    etoindex = Poly([-shift/gain, 1/gain])
    etoindexint = Int∘round∘etoindex
    indextoe, etoindex, etoindexint
end

function transformations(detector::ChannelDetector, nbins=5000)
    hist, edges = histogram(detector, nbins=nbins)
    transformations(edges)
end
