using .JSort

using SharedArrays

using JLD
using DelimitedFiles
using Interpolations
using Polynomials
import Base.show
import YAML

function constructzrange(path)
    data = readdlm(path, skipstart=2)
    itp = LinearInterpolation(convert(Array{Float64, 1}, 1000data[:, 1]),
                              convert(Array{Float64, 1}, data[:, 2]),  # Convert to keV
                              extrapolation_bc=Line())
    return itp
end

function readgainshift(path::Nothing)
    gainshift = [Float32.([0.0, 1.0]) for f in 1:8, b in 1:8]
end

function readgainshift(path)
    println(path)
    gainshift, header = readdlm(path, ',', header=true)
    # Check the file
    b = f = 1
    gainshiftarray = [Float32[] for f in 1:8, b in 1:8]
    for (i, row) in enumerate(eachrow(gainshift))
        # if Int(row[1]) ≠ f || Int(row[2]) ≠ b
        #     error("Gainshift file is corrupt at row $i = $row")
        # end
        f = Int(row[1])
        b = Int(row[2])
        gainshiftarray[f, b] = row[3:end]
    #     f += 1
    #     if f > 8
    #         f = 1
    #         b += 1
    #     end
    end
    gainshiftarray
end

readgainshift_gamma(paths::Vector{Nothing}) = []
readgainshift_gamma(paths::Nothing) = []
function readgainshift_gamma(paths::Vector{T}) where T
    [readgainshift_gamma(path) for path in paths]
end

function readgainshift_gamma(path::AbstractString)
    gainshift, header = readdlm(path, ',', header=true)
    gainshiftarray = zeros((33, 3))
    channels = Set(collect(1:33))
    for (i, row) in enumerate(eachrow(gainshift))
        try
            channel = round(Int, row[1])
            gainshiftarray[channel, :] = row[2:end]
            pop!(channels, channel)
        catch e
            println(e)
            error("Gamma gainshift file corrupt at row $i = $row")
        end
    end
    if length(channels) > 0
        @warn "Channels $channels are unspecified $path"
    end
    gainshiftarray
end

struct Parameters
    parameterpath::String
    raw  # The YAML file
    calibrator::Dict{Symbol, AbstractArray{Poly}}
    isotope
    datafiles
    exfromeΔe::Array{Float64, 2}
    excorr::Array{Float64, 2}
    thickrange::Array{Float64, 1}
    tnaicorresi::Array{T, 1} where T<:Number
    eΔerect
    percent
    range::Interpolations.Extrapolation{Float64,1,Interpolations.GriddedInterpolation{Float64,1,Float64,Interpolations.Gridded{Interpolations.Linear},Tuple{Array{Float64,1}}},Interpolations.Gridded{Interpolations.Linear},Interpolations.Line{Nothing}}
    cores
    savepath::String
    gates::Dict{Symbol,Gate}

    function Parameters(parameterpath::T) where T <: AbstractString
        join(path) = joinpath(splitdir(parameterpath)[1], path)
        join(path::Nothing) = nothing
        parameters    = YAML.load_file(parameterpath)
        isotope       = getparameter(parameters, "isotope")
        datafiles     = getparameter(parameters, "datafiles", T=Vector{String}) .|> join
        exfromeΔe     = getparameter(parameters, "ex from ede", repeat([0.0 0.0 0.0], 8))
        excorr        = getparameter(parameters, "excitation energy correction", repeat([0.0 1.0], 8))
        thickrange    = getparameter(parameters, "thick range", T=Array{Float64})
        tnaicorresi   = getparameter(parameters, "tnai corr enai", [0, 0, 1, 0])
        eΔerect       = getparameter(parameters, "ede rectangle")
        percent       = getparameter(parameters, "percent to read", 1.0)
        cores         = getparameter(parameters, "cpu cores", 1)
        savepath      = getparameter(parameters, "save path", "sirius") |> join

        # Calibration
        egainshift  = readgainshift(join(parameters["e gainshift file"]))
        Δegainshift = readgainshift(join(parameters["de gainshift file"]))
        γgainshift  = parameters["gamma gainshift file"] .|> join |> readgainshift_gamma
        tgainshift  = parameters["time gainshift file"]  .|> join |> readgainshift_gamma
        calibrator  = Dict{Symbol, Calibrator}(
            :e  => combinegainshift(egainshift, 1:8, 1:8),
            :Δe => combinegainshift(Δegainshift, 1:8, 1:8, initial=[0.0, 2.5]),
            :γ  => combinegainshift(γgainshift, 1:33),
            :t  => combinegainshift(tgainshift, 1:33, initial=[0.0, 1.0])
        )
        # Gates # Fix case for empty gates
        gates = getgates(parameters["gates"])
        range = getparameter(parameters, "range file") |> join |> constructzrange

        new(parameterpath, parameters, calibrator, isotope, datafiles, exfromeΔe, excorr, thickrange,
            tnaicorresi, eΔerect, percent, range, cores, savepath, gates)
    end
end

function show(io::IO, a::Parameters)
    println("Parameters for $(a.isotope)")
    println("Loaded from $(a.parameterpath).")
    println("e-calibration coefficients:")
    show(io, a.calibrator[:e])
    println("\nΔe-calibration coefficients:")
    show(io, a.calibrator[:Δe])
    println("\ngamma calibration coefficients:")
    show(io, a.calibrator[:γ])
    println("\ntime calibration coefficients:")
    show(io, a.calibrator[:t])
    println("Gates:")
    for gate in a.gates
        println("$(gate.first) : $(gate.second.low) to $(gate.second.high)")
    end
end

function show(io::IO, a::AbstractVector{T} where T<: Poly)
    for i in eachindex(a)
        print("$i: ")
        println(a[i])
    end
end

function show(io::IO, a::AbstractMatrix{T} where T<: Poly)
    for row in 1:size(a, 1)
        for col in 1:size(a, 2)
            print("$row, $col: ")
            println(a[row, col])
        end
    end
end


function getparameter(parameters, parameterkey, defaultvalue; constraint=nothing)
    if parameterkey ∉ keys(parameters)
        parameter = defaultvalue
    else
        parameter = parameters[parameterkey]
        # Reshape arrays
        if typeof(defaultvalue) <: Array{Float64}
            parameter = reshape(parameter, size(defaultvalue))
        else
            parameter = convert(typeof(defaultvalue), parameter)
        end
        if !isnothing(constraint)
            if !constraint(parameter)
                error("Parameter $parameterkey does not satisfy $constraint")
            end
        end
    end
    return parameter
end

function getparameter(parameters, parameterkey; T=nothing, constraint=nothing)
    if parameterkey ∉ keys(parameters)
        error("Parameter $parameterkey is required")
    else
        parameter = parameters[parameterkey]
        if !isnothing(T)
            parameter = convert(T, parameter)
        end
        if !isnothing(constraint)
            if !constraint(parameter)
                error("Parameter $parameterkey does not satisfy $constraint")
            end
        end
    end
    return parameter
end


function getpath(parameters::Parameters, key::Symbol)
    basepath = splitdir(parameters.parameterpath)[1]
    path = joinpath(basepath, getfield(parameters, key))
end

function getpath(parameters::Parameters, key)
    basepath = splitdir(parameters.parameterpath)[1]
    filepath = parameters.raw[key]
    path = joinpath(basepath, filepath)
end

### HANDLE GAINSHIFT ###
function combinegainshift(gainshifts::AbstractVector, channels;
                          initial = [0.0, 5.0],
                          docalibrate = true)::Vector{Poly}
    #@assert length(channels) == length(gainshift)
    # Combine all gainshift operations into a single
    # loop for a given channel. Gives a 3x speedup
    # Try to unroll the loop
    calibrators = Poly[]
    for channel in channels
        p = Poly(initial)
        if docalibrate
            for gainshift in gainshifts
                p = Poly(gainshift[channel, :])(p)
            end
        end
        # c = coeffs(p)
        # n = length(c)
        # f = let n=n, c=c
            # x -> begin
                # s = 0.0
                # @inbounds for i in 1:n
                    # s += c[i]*x^(i-1)
                # end
                # s
            # end
        # end
        push!(calibrators, p)
    end
    calibrators
end

function combinegainshift(gainshifts::AbstractMatrix, channelsf, channelsb;
                          initial = [0.0, 5.0],
                          docalibrate = true)::Matrix{Poly}
    combinegainshift([gainshifts], channelsf, channelsb, initial=initial)
end
function combinegainshift(gainshifts::AbstractVector, channelsf, channelsb;
                          initial = [0.0, 5.0],
                          docalibrate = true)::Matrix{Poly}
    # Combine all gainshift operations into a single
    # loop for a given channel. Gives a 3x speedup
    # Try to unroll the loop
    calibrators  = Matrix{Poly}(undef, length(channelsf), length(channelsb))
    for f in channelsf, b in channelsf
        p = Poly(initial)
        if docalibrate
            for gainshift in gainshifts
                p = Poly(gainshift[f, b])(p)
            end
        end
        # TODO why not just use Poly(x)??
        # c = coeffs(p)
        # n = length(c)
        # func = let n=n, c=c
            # x -> begin
                # s = 0.0
                # @inbounds for i in 1:n
                    # s += c[i]*x^(i-1)
                # end
                # s
            # end
        # end
        calibrators[f, b] = p
    end
    calibrators
end

function savegainshift_gamma(path, coefficients)
    println(coefficients)
end

function savegainshift_gamma(parameters::Parameters, coefficients)
    path = getpath(parameters, "gamma gainshift file")
    savegainshift_gamma(path, coefficients)
end

function savegainshift_gamma(parameters::Parameters, filename::String, coefficients)
    path = getpath(parameters, :savepath)
    path = joinpath(path, filename)
    # pad different lengthed polynomials
    N = maximum(map(length, coefficients))
    coefficients = map(x -> begin
            diff = length(x) - N
            if diff < 0
                x = [x..., zeros(abs(diff))...]
            else
                x
            end
        end, coefficients)

    open(path, "w") do io
        header = ["channel" "shift" "gain" "sqgain"]
        @show coefficients
        writedlm(io, [header, coefficients...], ',', header=true)
    end
end


function loadlabr(parameters::Parameters)::Vector{TrueEvent}
    return loadlabr(parameters.savepath)
end

function loadlabr(path::AbstractString)::Vector{TrueEvent}
    fhandle = open(joinpath(path, "labr.bin"))
    data = Base.read(fhandle) |> x -> reinterpret(UInt32, x)
    i = 1
    labr = TrueEvent[]
    while i ≤ length(data)
        f, b = data[i:i+1]
        i += 2
        de, e = convert.(Float32, data[i:i+1])
        i += 2
        labrcount = data[i]
        if labrcount + i >= length(data)
            println("Reading failed!")
        end
        i += 1
        adctdcs = ADCTDC[]
        for j in 1:labrcount
            channel = data[i]
            i += 1
            adc, tdc = convert.(Float32, data[i:i+1])
            push!(adctdcs, ADCTDC(channel, adc, tdc))
            i += 2
        end
        push!(labr, TrueEvent(f, b, de, e, adctdcs))
    end
    labr
end

function savecoefficients(parameters, path, coefficients, channels)
    all_coefficients = [[0.0, 1.0, 0.0] for i in 1:33]
    all_coefficients[channels] = coefficients
    # Tag each coefficient with corresponding channel
    all_coefficients = [[i, val...] for (i, val) in enumerate(all_coefficients)]
    savegainshift_gamma(parameters, path * ".csv", all_coefficients)
end

function savecoefficients_eΔe(parameters, coefficients_e, coefficients_Δe, path::Tuple{String, String}=("egainshift.csv", "degainshift.csv"))
    # Coefficients formated as  front, back, x0   x1   y0   y1
    savepath = getpath(parameters, :savepath)
    spath = joinpath(savepath, path[1])
    open(spath, "w") do io
        header = ["front" "back" "shift" "gain"]
        @show coefficients_e
        writedlm(io, [header, coefficients_e...], ',', header=true)
    end

    spath = joinpath(savepath, path[2])
    open(spath, "w") do io
        header = ["front" "back" "shift" "gain"]
        @show coefficients_Δe
        writedlm(io, [header, coefficients_Δe...], ',', header=true)
    end
end

function readede(path::String; T=Float32)
    data = open(path, "r") do io
        len = read(io, T) |> Int
        data = Matrix{T}(undef, len, 2)
        read!(io, data)
        data
    end
end
