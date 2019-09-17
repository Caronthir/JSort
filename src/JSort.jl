module JSort
export processfile, Event, sortfile, save, Parameters, OMatrix, calibrate, xtoindex, measurequality, increment!, creatematrices, save, loadlabr, makeeΔe, MAX_E, MAX_ΔE, NUMBINS_E, gainshift!, makelabr, linearalign, goodchannels, findpeaks, gainshift, leastsquares, localalign, splitregion, featurealign, boxbounds, arraybox, splitregions2d, featureoverlap2d, featurealign2d, indextox, indextoy, xrange, yrange, alignspectra, ucalibrate
include("event.jl")
include("filehandling.jl")
include("reader.jl")
include("matrices.jl")
include("sorter.jl")
include("manager.jl")
include("alignment.jl")
include("alignment2d.jl")
end
