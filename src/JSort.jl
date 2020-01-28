module JSort
export Event, Parameters, Reader,
       linearalign, goodchannels, findpeaks, gainshift,
       leastsquares, localalign, splitregion, featurealign, boxbounds,
       arraybox, splitregions2d, featureoverlap2d, featurealign2d,
       alignspectra, bin, binhist, aligntime,
       gate, Gate, histogram_edges, histogram, L1, L2, χ², Weight, ThresholdSum,
       ThresholdArea, featurealign2d_weighted, annotate_plot_2d, compare_2d,
    L1l, L2l, neighbourreduce, plotresiduals, compare_methods, splitmetric,
    ParticleDetector, ParticleCalibrator, read, readparticledetector,
    fill!, combine!, addcalibrator!

include("gate.jl")
include("calibrator.jl")
include("detectors.jl")
include("parameters.jl")
include("metrics.jl")
include("weights.jl")
include("regression.jl")
include("reader.jl")
include("manager.jl")
include("alignment.jl")
include("alignment2d.jl")
include("histogram.jl")
include("timealignment.jl")

#include("BetheBloch/BetheBloch.jl")
#using .BetheBloch
#export Material, Particle, Banana

end
