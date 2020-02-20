using .JSort
using DelimitedFiles: readdlm

function combinepeaks(calib::AbstractVector{V} where V<:ParticleCalibrator, fnames::AbstractVector{T} where T<:AbstractString)
    # We want to backtransform the peaks so that all peaks
    # are given in original channel numbers, which are the same for all
    # experimental data.
    calib = invert.(calib)
    # Load file
    peaks_exp = [Float64[] for f in 1:8, b in 1:8]
    peaks_kinz = [Float64[] for f in 1:8, b in 1:8]
    for (i, file) in readdlm.(fnames, ',') |> enumerate
        for line in eachrow(file)
            front, e, qkinz = line
            front = Int(front) + 1
            for back in 1:8
              push!(peaks_exp[front, back], calib[i][front, back](e))
              push!(peaks_kinz[front, back], qkinz)
            end
        end
    end

    newcalib = ParticleCalibrator()
    for front in 1:8, back in 1:8
        @assert length(peaks_exp[front, back]) > 1 "Provide at least 2 peaks"
        coeff = leastsquares(peaks_exp[front, back], peaks_kinz[front, back]; order=1)
        newcalib[front, back] = Poly(coeff, calib[1][1,1].var)
    end
    return newcalib
end
