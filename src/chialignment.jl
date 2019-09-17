#= Ideas that didn't really work in pratice
=#
function localalign(target::Vector{T}, reference::Vector{T}, roi;
                    numregions=10, method=:peak) where T
    @assert length(target) == length(reference)
    # TODO I must only shift a window around the region of interest
    # If not, then the outskirts will bias the minimization, and the local
    # inverse transformation will be off
    # - What did I mean by this??

    points = []
    for (start, stop) in splitregion(roi, numregions)
        println("Working on ", start, ", ", stop)
        region = first(roi) .- 1 .+ [start:stop;]
        reference_point, target_point = try regionalign(target, reference, region)
        catch e
            println("==================")
            println(e)
            if isa(e, BoundsError)
                continue
            end
            throw(e)
        end

        # Shift points back
        reference_point += start + first(roi)
        target_point += start + first(roi)
        push!(points, [reference_point, target_point])
    end
    global tracecounter
    tracecounter = tracecounter + 1
    points
end
function regionalign(target::Vector{T}, reference::Vector{T}, roi, method=:peak,
                     smoothing=[31, 4]) where T
    @assert method ∈ [:peak, :random]
    signal = pyimport("scipy.signal")
    optim = pyimport("scipy.optimize")
    smooth(x) = signal.savgol_filter(x, smoothing...)[roi]

    refpoint = 0
    ref = smooth(reference)
    if method == :peak
        refpoint = findpeaks(ref, numpeaks=1)[1]
    elseif method == :random
        refpoint = rand(roi)
    end

    # Minimize the reference against to target to obtain the transformation
    # R→T
    shift = minimize(target, reference, roi, smooth)
    #corr = signal.correlate(target[roi], reference[roi])
    #plt = plot(corr)
    #display(plt)
    #shift = - (length(roi) - argmax(corr))
    gain = 1.0
    targetpoint = refpoint + shift |> x -> round(Int, x)

    shifted = gainshift(reference, shift, gain)
    #tracealignment(reference, target, shifted, refpoint, targetpoint, smooth, roi)

    refpoint, targetpoint
end

function log10transform(V)
    G = V[1:end]
    G[G .> 0] .= log10.(G[G .> 0])
    G[G .≤ 0] .= NaN
    G
end

tracecounter = 0
function tracealignment(reference::Vector{T}, target::Vector{T}, shifted::Vector{T}, refpoint, targetpoint, smoother, roi) where T
    tcol = :blue
    rcol = :red
    scol = :green
    bcol = :purple
    ycol = :yellow

    ref = smoother(reference)
    smoothtarget = smoother(target)
    plt = plot(reference[roi], label="Reference", alpha=0.3, color=rcol, seriestype=:steppre)
    plt = plot!(ref, label="Reference, smooth", color=rcol)
    plt = scatter!([refpoint], [ref[refpoint]],  markershape=:dtriangle, markercolor=rcol,
                   markerstrokewidth=0, label="Peak, smooth")

    plt = plot!(target[roi], label="Target", alpha=0.3, color=tcol, seriestype=:steppre)
    plt = plot!(smoothtarget,  label="Target, smooth", color=tcol)

    plt = plot!(shifted[roi], label="Shifted reference", alpha=0.3, color=scol, seriestype=:steppre)
    plt = plot!(smoother(shifted), label="Shifted reference, smooth", color=scol)

    plt = scatter!([targetpoint], [smoothtarget[targetpoint]],  markershape=:dtriangle, markercolor=bcol,
                   markerstrokewidth=0, label="Shifted peak onto target")

    l = @layout [  a{0.15h}; b{0.15h}; c]
    G = reference[roi] |> x -> convert(Vector{Float64}, x) |> log10transform
    P = target[roi] |> x -> convert(Vector{Float64}, x) |> log10transform
    p1 = heatmap(G')
    p1 = plot!([refpoint], [1], markershape=:dtriangle, markercolor=rcol, markerstrokewidth=0)
    p2 = heatmap(P')
    p2 = plot!([targetpoint], [1], markershape=:dtriangle, markercolor=tcol, markerstrokewidth=0)
    p = plot(p1, p2, plt, layout=l, yticks=nothing, colorbar=false, margin=0.2mm, size=(1000, 1000),
            legend=false)
    #p.subplots[1].attr[:xaxis][:ticks] = nothing
    p.subplots[1].attr[:framestyle] = :none

    global tracecounter
    savefig(p, "logplots/strip_$(tracecounter)_region_$(first(roi)).png")
end
