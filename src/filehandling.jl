using .JSort

using SharedArrays
using JLD
using DelimitedFiles
using Interpolations
import YAML

function constructzrange(path)
    data = readdlm(path, skipstart=2)
    itp = LinearInterpolation(convert(Array{Float64, 1}, 1000data[:, 1]),
                              convert(Array{Float64, 1}, data[:, 2]),  # Convert to keV
                              extrapolation_bc=Line())
    return itp
end

function readgainshift(path)
    gainshift, header = readdlm(path, ',', header=true)
    # Check the file
    b = f = 1
    gainshiftarray = zeros((8, 8, 3))
    for (i, row) in enumerate(eachrow(gainshift))
        if row[1] ≠ b || row[2] ≠ f
            error("Gainshift file is corrupt at row $i = $row")
        end
        gainshiftarray[b, f, :] = row[3:end]
        f += 1
        if f > 8
            f = 1
            b += 1
        end
    end
    gainshiftarray
end

struct Parameters
    gainshift::Dict{Symbol, Array{Float64, 3}}
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
    function Parameters(parameterpath::T) where T <: AbstractString
        parameters    = YAML.load_file(parameterpath)
        isotope       = getparameter(parameters, "isotope")
        datafiles     = getparameter(parameters, "datafiles", T=Array{String})
        exfromeΔe     = getparameter(parameters, "ex from ede", repeat([0.0 0.0 0.0], 8))
        excorr        = getparameter(parameters, "excitation energy correction", repeat([0.0 1.0], 8))
        thickrange    = getparameter(parameters, "thick range", T=Array{Float64})
        tnaicorresi   = getparameter(parameters, "tnai corr enai", [0, 0, 1, 0])
        eΔerect       = getparameter(parameters, "ede rectangle")
        percent       = getparameter(parameters, "percent to read", 1.0)
        cores         = getparameter(parameters, "cpu cores", 1)
        savepath      = getparameter(parameters, "save path", "sirius")
        gainshift = Dict(:e  => readgainshift(parameters["e gainshift file"]),
                         :Δe => readgainshift(parameters["de gainshift file"]))
        range = constructzrange(getparameter(parameters, "range file"))
        new(gainshift, isotope, datafiles, exfromeΔe, excorr, thickrange,
            tnaicorresi, eΔerect, percent, range, cores, savepath)
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

function save(labrs::Array{LaBrEvent}, path::T) where T<:AbstractString
    # Terrible performance. Concoct some smarter way to save data
    !isdir(path) && mkdir(path)
    fhandle = open(joinpath(path, "labr.bin"), "w")
    for labr in labrs
        header, body = serialize(labr)
        write(fhandle, header) 
        for b in body
            write(fhandle, b...)
        end
    end
end

function loadlabr(path)::Array{TrueEvent, 1}
    fhandle = open(joinpath(path, "labr.bin"))
    data = read(fhandle) |> x -> reinterpret(UInt32, x)
    i = 1
    labr = TrueEvent[]
    while i <= length(data)
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
