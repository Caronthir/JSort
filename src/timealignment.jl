using .JSort
using PyCall
import PyPlot; plt = PyPlot;
using Polynomials

#TODO Can maybe improved by fitting the actual exponential decay
function aligntime(spectra, edges; referenceindex=1,
                   prominence=50, distance=290, plot=false,
                   prompt=400, wavelength=45)
    signal = pyimport("scipy.signal")
    coefficients = Vector{Float64}[]
    aligned = Vector{eltype(spectra[1])}[]
    allpeaks = Vector{Float64}[]

    shift = first(edges)
    gain = step(edges)
    indextot = Poly([shift, gain])
    ttoindex = Poly([-shift/gain, 1/gain])
    distance = ttoindex(distance)

    for spectrum in spectra
        # Find the main peak
        peak = findpeaks(spectrum, numpeaks=1)[1]

        # Search for all smaller peaks to the right for the main peak
        region = (peak+40):length(spectrum)
        peaks, prop = signal.find_peaks(spectrum[region],
                                        prominence=prominence,
                                        distance=distance)
        peaks .+= first(region) # shift back
        peaks = [peak, peaks...]
        push!(allpeaks, peaks)
    end

    for (chn, peaks) in allpeaks |> enumerate
        # Last peak is often ugly
        N = length(peaks) - 1
        N = if N > 10 10 else N end
        peaks = peaks[1:N]
        peaks .= indextot.(peaks)
        # Main peak at 400 ms, next follows with distance 45 ms
        fact = [prompt, [prompt + wavelength*i for i in 1:N-1]...]
        @show fact
        @show peaks
        coeffs = leastsquares(peaks, fact)
        push!(coefficients, coeffs)
        if plot
            plt.scatter(peaks, repeat([chn+0.5], length(peaks)), color="r", marker="v")
        end
    end

    # flag = false
    # for (i, (peaks, spectrum)) in zip(allpeaks, spectra) |> enumerate
    #     if i == referenceindex
    #         push!(coefficients, [0.0, 1.0])
    #         push!(aligned, spectrum)
    #         continue
    #     end
    #     # Use as many peaks as possible for alignment
    #     refpeaks = allpeaks[referenceindex]
    #     N = min(length(peaks), length(refpeaks))
    #     # The last peak is often ugly. Remove it
    #     N = N-1
    #     refpeaks = refpeaks[1:N]
    #     peaks = peaks[1:N]
    #     coeffs = leastsquares(peaks, refpeaks)
    #     #push!(aligned, gainshift(spectrum, coeffs...))

    #     # Convert to time basis
    #     refpeaks .= indextot.(refpeaks)
    #     peaks .= indextot.(peaks)
    #     coeffs = leastsquares(peaks, refpeaks)
    #     push!(coefficients, coeffs)
    #     if plot
    #         plt.scatter(peaks, repeat([i+0.5], length(peaks)), color="r", marker="v")
    #         if !flag
    #             plt.scatter(refpeaks, repeat([referenceindex+0.5], length(refpeaks)), color="r", marker="v")
    #             flag = true
    #         end
    #     end
    # end
    coefficients
end

#=
# t: Uncorrected time
# e: Calibrated energy [keV]
# e_si: Calibrated SiRi back energy [keV]
=#
function chirp(time, e, e_si)
    corr[0]
end
