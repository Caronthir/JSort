using .JSort
import PyPlot; const plt = PyPlot
import PyCall
#using Interactive

function GatePicker(matrix, x, y)
    pushfirst!(PyVector(pyimport("sys")."path"), "/home/erdos/gits/JSort/src/")
    interactivegating = pyimport("interactivegating")
    GP = interactivegating.GatePicker
    return GP(matrix, x, y)
end

function GatePicker(det::ParticleDetector; nbins=1000)
    (matrix, (x, y)) = histogram(det, nbins=nbins)
    gp = GatePicker(matrix, x, y)
    return () -> Gate2D(gp.extent...)
end

@enum Action Continue Quit Align
struct ActionEvent
    action::Action
    data::Vector{Any}
end

function ActionEvent(input::AbstractString)
    if input == "q"
        return ActionEvent(Quit, [])
    end
    tokens = split(input, r"\s+")
    if length(tokens) ≥ 2
        try
            data = parse.(Float64, tokens[1:2])
            for token in tokens[3:end]
                key, val = split(token, r"\s*=\s*")
                push!(data, Symbol(key) => parse(Float64, val))
            end
            return ActionEvent(Align, data)
        catch
            return ActionEvent(Continue, [])
        end
    end
    return ActionEvent(Continue, [])
end

quit(event::ActionEvent) = event.action == Quit
align(event::ActionEvent) = event.action == Align

function pickpgate(det::ParticleDetector, parameters::Parameters; nbins=1000)
    quick = QuickReader(parameters)
    read(quick)
    ungated = ChannelDetector(deepcopy(quick))
    char = ""
    default = GammaCalibrator([-913.5, 5.6])
    default(ungated)
    while true
        g = GatePicker(det, nbins=nbins)
        # Pause to let Python write to g
        event = ActionEvent(readline())
        aligner = (x...) -> ()
        if quit(event)
            break
        elseif align(event)
            aligner = (gated) -> begin
                calib = ialign(ungated, event.data[1]:event.data[2], event.data[3:end]...)
                calib(gated)
            end
        end
        @show event

        read(quick, g())
        det2 = ChannelDetector(quick)
        default(det2)
        aligner(det2)
        det2
        plot(det2)
        plot(det2 |> sum)
        plot(det2[1])
    end
end

function ialign(detector::ChannelDetector, region; nbins=5000, kwargs...)
    full, edges = histogram(detector, nbins=nbins)
    chns = goodchannels(detector)

    signal = pyimport("scipy.signal")
    smoother(x) = signal.savgol_filter(x, 51, 4)

    coeff, _ = alignspectra(full, edges, lowregion=region, lowsmoother=smoother,
                            lownumregions=get(kwargs, :nregions, 2),
                            lowregionwidth=get(kwargs, :swidth, 200),
                            lowsearchwidth=get(kwargs, :lwidth, 60),
                            plot=get(kwargs, :plot, false), order=1)
    coefficients = repeat([0.0 1.0], outer=(33, 1))
    coefficients[chns, :] .= hcat(coeff...)'
    calib = GammaCalibrator(coefficients)
end
