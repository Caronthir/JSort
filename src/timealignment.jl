using .JSort
using PyCall

#TODO Can maybe improved by fitting the actual exponential decay
function aligntime(spectra, edges; referenceindex=1,
                    prominence=50, distance=15)
    signal = pyimport("scipy.signal")
    coefficients = Vector{Float64}[]
    aligned = Vector{eltype(spectra[1])}[]
    allpeaks = Vector{Int}[]
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

    for (i, (peaks, spectrum)) in zip(allpeaks, spectra) |> enumerate
        if i == referenceindex
            push!(coefficients, [0.0, 1.0, 0.0])
            push!(aligned, spectrum)
            continue
        end
        # Use as many peaks as possible for alignment
        refpeaks = allpeaks[referenceindex]
        N = min(length(peaks), length(refpeaks))
        # The last peak is often ugly. Remove it
        N = N-1
        refpeaks = refpeaks[1:N]
        peaks = peaks[1:N]
        coeffs = leastsquares(peaks, refpeaks)
        push!(aligned, gainshift(spectrum, coeffs...))

        # Convert to time basis
        shift = first(edges)
        gain = step(edges)
        refpeaks = refpeaks .* gain .+ shift
        peaks = peaks .* gain .+ shift
        coeffs = leastsquares(peaks, refpeaks)
        push!(coefficients, coeffs)
    end
    aligned, coefficients
end

#=
# t: Uncorrected time
# e: Calibrated energy [keV]
# e_si: Calibrated SiRi back energy [keV]
=#
function chirp(time, e, e_si)
    corr[0]
end
