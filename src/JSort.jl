module JSort
export Parameters,
       linearalign, goodchannels, findpeaks, gainshift,
       leastsquares, localalign, splitregion, featurealign, boxbounds,
       arraybox, splitregions2d, featureoverlap2d, featurealign2d,
       alignspectra, bin, binhist, aligntime,
       histogram_edges, histogram, L1, L2, χ², Weight, ThresholdSum,
       ThresholdArea, featurealign2d_weighted, annotate_plot_2d, compare_2d,
    L1l, L2l, neighbourreduce, plotresiduals, compare_methods, splitmetric,
    ParticleDetector, ParticleCalibrator, read, readparticledetector,
    fill!, combine!, addcalibrator!, CBuffer, CBuffer2D, read, write!, next!, close,
    open, isgood, plot, calibrateall, invert, getindex

include("gainshift.jl")
export gainshift

include("./BetheBloch/BetheBloch.jl")
export Banana

include("./buffer/buffer.jl")
include("gate.jl")
export Gate, center, gate, center, Point, Gate2D, AbstractGate

include("cascade.jl")
export GammaCascade, extend!

include("levels.jl")
export fetchlevels

include("calibrator.jl")
export TimeFunction, TimeCalibrator, GammaCalibrator
       TimeCorrector, AbstractCalibrator, AbstractCalibrator1D

include("peakcombine.jl")
export combinepeaks

include("detectors.jl")
export GammaDetector, TimeDetector, readres, Detector, ChannelDetector,
       align, align!

include("parameters.jl")
include("metrics.jl")

include("weights.jl")
export ThresholdSumShift, WellSum

include("regression.jl")
include("reader.jl")
export AbstractReader, Reader, Event

include("quickreader.jl")
export QuickReader

include("alignment.jl")

include("alignment2d.jl")
export automap

include("regionalign.jl")
export regionpeaks, regionalign, backshift, plothist

include("iterativealign.jl")
export iterativealign, GammaCascade

include("histogram.jl")
include("timealignment.jl")

include("interactive.jl")
export GatePicker, pickpgate


#include("BetheBloch/BetheBloch.jl")
#using .BetheBloch
#export Material, Particle, Banana

end
