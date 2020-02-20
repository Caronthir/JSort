using .JSort
import Base: read, read!, show, length, sum
using DataFrames
import PyPlot; const plt = PyPlot
using LaTeXStrings

macro tostring(arg)
    string(arg)
end

abstract type AbstractDetector end
abstract type AbstractDetector1D <: AbstractDetector end
abstract type AbstractDetector2D <: AbstractDetector end

length(d::AbstractDetector1D) = length(d.ref[])

function read(path::AbstractString, V::Type{W}; T=UInt32) where W <: AbstractDetector
    detector = open(path, "r") do io
        read(io, V, T=T)
    end
end

function read(io::IO, V::Type{<:AbstractDetector1D}; T=UInt32)
    len = stat(io).size / sizeof(T) |> Int
    data = Vector{Float32}(undef, len)
    for i in 1:len
        data[i] = read(io, T) |> Float32
    end
    pattern = r".*chn(?<c>\d+).*"
    m = match(pattern, io.name)
    chn = if isnothing(m) 0 else parse(Int, m[:c]) end
    V(data, chn)
end

function read(io::IO, V::Type{W} where W <: AbstractDetector2D; T=UInt32)
    # Assumes binary format int UInt32
    len = stat(io).size / 2sizeof(T) |> Int
    data = Matrix{Float32}(undef, len, 2)
    for i in 1:len
        data[i, 1] = read(io, T) |> Float32
        data[i, 2] = read(io, T) |> Float32
    end

    pattern = r".*b(?<b>\d)f(?<f>\d).*"
    m = match(pattern, io.name)
    if isnothing(m)
        pattern = r".*f(?<f>\d).*"
        m = match(pattern, io.name)
        if isnothing(m)
            f = b = 0
        else
            f = parse(Int, m[:f])
            b = 0
        end
    else
        f, b = parse.(Int, (m[:f], m[:b]))
    end
    V(data[:, 1], data[:, 2], f, b)
end

function write(path::AbstractString, V::AbstractDetector)
    open(path, "w") do io
        write(io, V)
    end
end

function write(io::IO, V::AbstractDetector1D)
    write(io, V.ref[])
end

function write(io::IO, V::AbstractDetector2D)
    for i in eachindex(V.Xref[])
        write(io, V.Xref[][i], V.Yref[][i])
    end
end

function histogram(X::AbstractDetector2D; nbins=1000)
    edges = histogram_edges(X.Xref[], X.Yref[], nbins=nbins)
    histogram(X, edges...), edges
end

function histogram(X::AbstractDetector2D, xedges, yedges)
    #histogram(X.Xref[], X.Yref[], xedges, yedges)
    histogram(X.Xref[] .+ rand(length(X.Xref[])),
              X.Yref[] .+ rand(length(X.Yref[])),
              xedges, yedges)
end

function histogram(X::AbstractDetector2D, edges)
    histogram(X.Xref[], X.Yref[], edges...)
end

function histogram(X::AbstractDetector1D, edges::AbstractRange)
    # Why do I need to add a random number to remove histogram artifacts?
    histogram(X.ref[] .+ rand(length(X.ref[])), edges)
    #histogram(X.ref[], edges)
end
function histogram(X::AbstractDetector1D; nbins=1000)
    edges = range(minimum(X.ref[]), stop=maximum(X.ref[]), length=nbins)
    histogram(X, edges), edges
end

function (calib::AbstractCalibrator1D)(d::T; inplace::Bool=true)::T where {T<:AbstractDetector1D}
    poly = calib.poly[d.channel]
    if inplace
        d.ref[] .= poly.(d.ref[])
        return d
    else
        d = deepcopy(d)
        d.ref[] .= poly.(d.ref[])
        return d
    end
end

function show(io::IO, V::AbstractDetector1D)
    println("Channel: $(V.channel)")
    show(io, DataFrame(e=V.ref[]), allrows=false)
end
## Concrete implementations

struct ParticleDetector <: AbstractDetector2D
    e::Vector{Float32}
    Δe::Vector{Float32}
    front::Int
    back::Int
    Xref::Ref{Vector{Float32}}
    Yref::Ref{Vector{Float32}}
    function ParticleDetector(a::AbstractVector, b::AbstractVector, c::Number, d::Number)
        new(a, b, c, d, Ref(a), Ref(b))
    end
end
ParticleDetector(path::AbstractString; T=UInt32) = read(path, ParticleDetector, T=T)

# function write(io::IO, d::ParticleDetector)
#     for i in 1:length(d.e)
#         write(io, d.e[i], d.Δe[i])
#     end
# end

function show(io::IO, detector::AbstractDetector2D)
    println("Front: $(detector.front)\nBack: $(detector.back)")
    show(DataFrame(e=detector.Xref[], Δe=detector.Yref[]), allrows=false)
end

function (calib::ParticleCalibrator)(d::ParticleDetector)
    poly = calib.poly[d.front, d.back]
    if poly.var == :e
        map!(poly, d.e, d.e)
    elseif poly.var == :Δe
        map!(poly, d.Δe, d.Δe)
    else
        throw("Unable to deduce how to apply calibration. Specify `e` or `Δe`")
    end
end

function (calib::ParticleCalibrator)(d::ParticleDetector, var::Symbol)
    poly = calib.poly[d.front, d.back]
    if var == :e
        map!(poly, d.e, d.e)
    elseif var == :Δe
        map!(poly, d.Δe, d.Δe)
    else
        throw("Variable must be `e` or `Δe`, not $var")
    end
end

function calibrateall(path::AbstractString; epath="coefficients_e.txt", Δepath="coefficients_de.txt")
    j(x) = joinpath(path, x)
    ce = ParticleCalibrator(epath, var=:e)
    cde = ParticleCalibrator(Δepath, var=:Δe)

    for f in 1:8
        totale = Float32[]
        totalde = Float32[]
        for b in 1:8
            detector = ParticleDetector(j("edeb$(b)f$(f).bin"))
            ce(detector)
            cde(detector)
            append!(totale, detector.e)
            append!(totalde, detector.Δe)
        end
        write(j("edef$(f).bin"), ParticleDetector(totale, totalde, f, 0))
    end
end

function plot(d::ParticleDetector; ax=nothing)
    (fig, ax) = if isnothing(ax)
        plt.subplots()
    else
        ax.figure, ax
    end
    mat, edges = histogram(d, nbins=1000)
    #ax.scatter(d.e, d.Δe, marker=",", s=1, alpha=0.5)
    ax.pcolormesh(edges[1], edges[2], mat |> transpose .|> log10)
    ax.set_xlabel(L"E\quad [keV?]")
    ax.set_ylabel(L"\Delta E\quad [keV?]")
    ax.set_title("f$(d.front)b$(d.back)")
end

### Gamma Detector

struct GammaDetector <: AbstractDetector1D
    e::Vector{Float32}
    channel::Int
    ref::Ref{Vector{Float32}}
    GammaDetector(e::AbstractVector) = new(e, 0, Ref(e))
    GammaDetector(e::AbstractVector, chn::Int) = new(e, chn, Ref(e))
end
GammaDetector(path::AbstractString; T=Float32) = read(path, GammaDetector, T=T)
GammaDetector() = GammaDetector(Float32[])

function plot(D::AbstractDetector1D; ax=nothing, nbins=1000)
    (fig, ax) = if isnothing(ax)
        plt.subplots()
    else
        ax.figure, ax
    end
    hist, edges = histogram(D, nbins=nbins)
    ax.step(edges, hist)
    ax.set_yscale("log")
    ax.set_title("Channel $(D.channel)")
end

### Time 'Detector'

struct TimeDetector <: AbstractDetector1D
    t::Vector{Float32}
    channel::Int
    ref::Ref{Vector{Float32}}
    TimeDetector(t::AbstractVector) = new(t, 0, Ref(t))
    TimeDetector(t::AbstractVector, chn::Int) = new(t, chn, Ref(t))
end
TimeDetector(path::AbstractString; T=Float32) = read(path, TimeDetector, T=T)
TimeDetector() = TimeDetector(Float32[])

### Channel 'Detector'
struct ChannelDetector{T} <: AbstractDetector
    data::Vector{T}
end

function ChannelDetector(name::AbstractString, path="."; T=Float32)
    pattern = Regex("^" * name * "_chn(?<chn>\\d+).*")
    ChannelDetector(pattern, path, T=T)
end

function ChannelDetector(pattern::Regex, path="."; T=Float32)
    data = nothing
    for fname in readdir(path)
        m = match(pattern, fname)
        if !isnothing(m)
            name = splitext(basename(fname))[1]
            types = parsetype(name)
            Tcalib = types[1]
            if isnothing(data)
                data = Tcalib[Tcalib() for i in 1:33]
            end
            chn = parse(Int, m[:chn])
            data[chn] = Tcalib(joinpath(path, fname), T=T)
        end
    end
    ChannelDetector(data)
end

function show(io::IO, d::ChannelDetector)
    for chn in eachindex(d.data)
        println(io, "Channel $chn:\t$(length(d.data[chn]))")
    end
end

function plot(D::ChannelDetector; ax=nothing, nbins=1000)
    (fig, ax) = if isnothing(ax)
        plt.subplots()
    else
        ax.figure, ax
    end
    hists = Vector{Int}[]
    edges = histogram_edges(D.data[1].ref[], nbins=nbins)
    chns = goodchannels(D)
    for chn in chns
        hist = histogram(D.data[chn], edges)
        push!(hists, hist)
    end
    hists = cat(hists'..., dims=1)
    ax.pcolormesh(edges, 1:length(chns)+1, hists .|> log10)
    #ax.matshow(hists .|> log10)
    #ax.pcolormesh(1:length(edges), 1:length(chns), hists .|> log10)
    ax.set_yticks(1.5:length(chns)+1-0.5)
    ax.set_yticklabels(chns)
    ax
end

function goodchannels(D::ChannelDetector; deviationlimit=10)
    good_lengths = [chn for (chn, d) in enumerate(D.data) if length(d) > 100 && chn != 32 && chn != 30]
    lengths = length.(D.data[chn] for chn in good_lengths)
    means = mean.(D.data[chn].ref[] for chn in good_lengths)
    M = mean(means)
    deviation = @. 1/(lengths - 1)*√((means - M)^2)
    good_lengths[deviation .< deviationlimit]
end

function histogram(D::ChannelDetector; nbins=1000)
    chns = goodchannels(D)
    edges = histogram_edges(D.data[1].ref[], nbins=nbins)
    [histogram(D.data[chn], edges) for chn in chns], edges
end

function histogram(D::ChannelDetector, edges::AbstractRange)
    [histogram(D.data[chn], edges) for chn in goodchannels(D)]
end

function align(detector::ChannelDetector; nbins=1000, kwargs...)
    d = deepcopy(detector)
    calib, invcalib = align!(d, nbins=nbins; kwargs...)
    d, calib, invcalib
end

function align!(detector::ChannelDetector{GammaDetector}; nbins=1000, kwargs...)
    signal = pyimport("scipy.signal")
    smoother(x) = signal.savgol_filter(x, 51, 4)

    chns = goodchannels(detector)
    hist, edges = histogram(detector, nbins=nbins)
    coeff, invcoeff = alignspectra(hist, edges; highsmoother=smoother, kwargs...)
    def = if (length(coeff[1]) == 2)
        [0.0 1.0] else [0.0 1.0 0.0] end
    coefficients = repeat(def, outer=(33, 1))
    invcoefficients = copy(coefficients)

    coefficients[chns, :] .= hcat(coeff...)'
    calib = GammaCalibrator(coefficients)
    calib(detector)
    invcalib = GammaCalibrator(coefficients)
    calib, invcalib
end

function align!(detector::ChannelDetector{TimeDetector}; nbins=1000, kwargs...)
    chns = goodchannels(detector)
    hist, edges = histogram(detector, nbins=nbins)
    coeff = aligntime(hist, edges; kwargs...)
    coefficients = repeat([0.0 1.0], outer=(33, 1))
    coefficients[chns, :] .= hcat(coeff...)'

    calib = TimeCalibrator(coefficients)
    calib(detector)
    calib, calib
end

function (c::AbstractCalibrator1D)(d::ChannelDetector; inplace::Bool=true)::ChannelDetector
    if inplace
        map(c, d.data)
        return d
    else
        d = deepcopy(d)
        map(c, d.data)
        return d
    end
end

function sum(d::ChannelDetector{T}) where T
    T(cat([det.ref[] for det in d.data]...; dims=1))
end

getindex(X::ChannelDetector, i) = getindex(X.data, i)
getindex(X::ChannelDetector, I...) = getindex(X.data, I...)
setindex!(X::ChannelDetector, i) = setindex!(X.data, i)
setindex!(X::ChannelDetector, I...) = setindex!(X.data, I...)

### General loading
struct Detector <: AbstractDetector
    data::Array{Float32}
    meta::Dict{Symbol, Int}
end
function Detector(path::AbstractString; T=UInt32)
    data, meta = readres(path, T=T)
    meta = Dict(i => j for (i, j) in meta if !isnothing(j))
    Detector(data, meta)
end

function plot(D::Detector; ax=nothing, nbins=1000)
    if ndims(D.data) == 1
        plot1D(D, ax=ax, nbins=nbins)
    elseif ndims(D.data) == 2
        plot2D(D, ax=ax, nbins=nbins)
    end
end

function plot1D(D::Detector; ax=nothing, nbins=1000)
    (fig, ax) = if isnothing(ax)
        plt.subplots()
    else
        ax.figure, ax
    end
    hist, edges = histogram(D.data, nbins=nbins)
    ax.step(edges, hist)
    ax.set_title(maketitle(D.meta))
end

function plot2D(D::Detector; ax=nothing, nbins=1000)
    (fig, ax) = if isnothing(ax)
        plt.subplots()
    else
        ax.figure, ax
    end
    mat, edges = histogram(D.data, nbins=nbins)
    ax.pcolormesh(edges[1], edges[2], mat |> transpose .|> log10)
    ax.set_title(maketitle(D.meta))
end

function readres(path::AbstractString; T=UInt32)
    open(path, "r") do io
        ndim = filedimensions(io.name)
        data = if ndim == 1
            readres1D(io, T=T)
        elseif ndim == 2
            readres2D(io, T=T)
        else
            error("Impossible dimension $ndim")
        end
        data, parsefilename(path)
    end
end


function readres1D(io::IO; T=UInt32)::Vector{Float32}
    len = stat(io).size / sizeof(T) |> Int
    data = Vector{Float32}(undef, len)
    for i in 1:len
        data[i] = read(io, T) |> Float32
    end
    data
end

function readres2D(io::IO; T=UInt32)::Matrix{Float32}
    len = stat(io).size / 2sizeof(T) |> Int
    data = Matrix{Float32}(undef, len, 2)
    for i in 1:len
        data[i, 1] = read(io, T) |> Float32
        data[i, 2] = read(io, T) |> Float32
    end
    data
end

function parsefilename(fname::AbstractString)
    fb_pattern = r".*b(?<b>\d)f(?<f>\d).*"
    chn_pattern = r".*chn(?<chn>\d+).*"
    f_pattern = r".*f(?<f>\d).*"
    m = match(fb_pattern, fname)
    f = b = chn = nothing
    if isnothing(m)
        m = match(f_pattern, fname)
        if !isnothing(m)
            f = parse(Int, m[:f])
        end
    else
        f, b = parse.(Int, (m[:f], m[:b]))
    end
    m = match(chn_pattern, fname)
    if !isnothing(m)
        chn = parse(Int, m[:chn])
    end
    (:f => f, :b => b, :chn => chn)
end

function filedimensions(fname::AbstractString)::Int
    if occursin("edeb", fname)
        return 2
    end
    underscores = count("_", fname)
    all = if occursin(r"\d", fname) 1 else 0 end
    underscores - all + 1
end

function maketitle(dict::Dict{Symbol, Int})
    @show dict
    title = ""
    k = keys(dict)
    if :f in k
        title *= "f$(dict[:f])"
    end
    if :b in k
        title *= " b$(dict[:b])"
    end
    if :chn in k
        title *= "Channel $(dict[:chn])"
    end
    title
end

function parsetype(fname::AbstractString)::Vector{DataType}
    tokens = split(fname, "_")
    types = DataType[]
    for token in tokens
        if token == "eg"
            push!(types, GammaDetector)
        elseif token == "t" || token == "tc"
            push!(types, TimeDetector)
        end
    end
    types
end
