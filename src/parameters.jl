using .JSort
import YAML

function yamlvalue(yaml::Dict{String, Any}, key::String, default)
    if key ∉ keys(yaml)
        default
    else
        param = yaml[key]
    end
end

struct Parameters
    root::String
    readpath::String
    savepath::String
    datafiles::Vector{String}
    isotope::String
    calibrationpaths::Dict{Symbol, String}
    cores::Int
    gates::Dict{Symbol, Gate}
    calibrator::Dict{Symbol, AbstractCalibrator}
    percent::Float64
end

function Parameters(path::AbstractString)
    parameters::Dict{String,Any} = YAML.load_file(path)
    readkeys = Set{String}()
    function get(key::String, default)
        push!(readkeys, key)
        yamlvalue(parameters, key, default)
    end
    root     = splitdir(path)[1]
    savepath = get("save path", root)
    if startswith(savepath, "~")
        savepath = expanduser(savepath)
    else
        savepath = if savepath == root savepath else joinpath(root, savepath) end
    end

    readpath = get("read path", savepath)
    if startswith(readpath, "~")
        readpath = expanduser(readpath)
    else
        readpath = if readpath == savepath readpath else joinpath(root, readpath) end
    end
    isotope  = get("isotope", "")

    cores = get("cpu cores", 1)

    # Different calibration coefficients
    calibrationpaths = Dict{Symbol, String}(
        :e  => get("coefficients e",     "ecoefficients.csv"),
        :Δe => get("coefficients de",    "decoefficients.csv"),
        :γ  => get("coefficients gamma", "gcoefficients.csv"),
        :t  => get("coefficients time",  "tcoefficients.csv")
    )
    # Try to read from the coefficient paths, defaulting to identity polynomials
    calibrator = Dict{Symbol, AbstractCalibrator}(
        :e  => makepcalibrator(calibrationpaths[:e], :e),
        :Δe => makepcalibrator(calibrationpaths[:Δe], :Δe),
    )

    gates = getgates(get("gates", nothing))

    # Read the datafile paths
    datafiles = get("data files", nothing)
    if isnothing(datafiles)
        @warn "No data files will be processed!"
    else
        map!(x -> joinpath(readpath, x), datafiles, datafiles)
    end

    percent = get("percent to read", 100)

    for key in keys(parameters)
        if key ∉ readkeys
            @warn "The parameter $key was supplied but not used."
        end
    end

    Parameters(path, readpath, savepath, datafiles, isotope,
               calibrationpaths, cores, gates, calibrator,
               percent)
end

function makepcalibrator(path, var)::ParticleCalibrator
    if isfile(path)
        read(path, ParticleCalibrator, var=var)
    else
        ParticleCalibrator(var=var)
    end
end

### PARSE GATES ###
getgates(::Nothing) = Dict(:e => Gate(), :t => Gate())
function getgates(unparsedgates::AbstractDict)
    # Gates on the excitation energy
    ukeys = lowercase.(keys(unparsedgates))
    ekey =  ukeys ∩ ["ex", "excitation", "ex_int"]
    tkey = ukeys ∩ ["t", "time", "na_t_c"]

    gates = Dict{Symbol,Gate}(:ex => Gate(), :t => Gate())

    if length(ekey) == 1
        for gate in unparsedgates[ekey[1]]
            gates[:ex] = gates[:ex] ∩ parsegate(gate)
        end
    end

    if length(tkey) == 1
        gates[:t] = unparsedgates[tkey[1]] |> parsegate
    end

    if length(ukeys) > 0 && length(ekey ∪ tkey) == 0
        throw(IOError("Gates $ukeys not supported."))
    end
    gates
end

function parsegate(gate::AbstractDict)::Gate
    ukeys = gate |> keys .|> lowercase
    low = ukeys ∩ ["low", "l"]
    high = ukeys ∩ ["high", "h"]
    wrong = setdiff(ukeys, low ∪ high)
    if length(wrong) > 0
        throw(ArgumentError("Unsupported gate bounds: $wrong"))
    end

    if length(low) > 0 && length(high) > 0
        return Gate(low=gate[low[1]], high=gate[high[1]])
    elseif length(low) > 0
        return Gate(low=gate[low[1]])
    else
        return Gate(high=gate[high[1]])
    end
end

function parsegate(gate::AbstractArray)::Gate
    if length(gate) ≠ 2
        throw(ArgumentError("Prove an upper and lower limit for gate"))
    end
    return Gate(low=gate[1], high=gate[2])
end
