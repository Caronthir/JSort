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
    #datafiles::Vector{String}
    isotope::String
    gates::Dict{Symbol, Gate}
    calibrator::Dict{Symbol, AbstractCalibrator}
    percent::Float64
    badchannels::Set{Int}
    qkinz::Vector{Poly{Float32}}
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

    # Different calibration coefficients
    calibrationpaths = Dict{Symbol, Union{String, Vector{String}}}(
        :e  => get("coefficients e",     "coefficients_e.csv"),
        :Δe => get("coefficients de",    "coefficients_de.csv"),
        :γ  => get("coefficients gamma", "coefficients_g.csv"),
#        :γi => get("coefficients gamma inverse", "coefficients_g_i.csv"),
        :t  => get("coefficients time",  "coefficients_t.csv"),
        :tc => get("coefficients time correction", "coefficients_t_c.csv")
    )
    # Try to read from the coefficient paths, defaulting to identity polynomials
    calibrator = Dict{Symbol, AbstractCalibrator}(
        :e  => make(ParticleCalibrator, calibrationpaths[:e],  :e),
        :Δe => make(ParticleCalibrator, calibrationpaths[:Δe], :Δe),
        :γ  => make(GammaCalibrator,    calibrationpaths[:γ],  :eᵧ),
#        :γi => make(GammaCalibrator,    calibrationpaths[:γi], :eᵧ),
        :t  => make(TimeCalibrator,     calibrationpaths[:t],  :t),
        :tc => make(TimeCorrector,      calibrationpaths[:tc], :t)
    )

    gates = getgates(get("gates", nothing))

    # Read the datafile paths
    # datafiles = get("data files", nothing)
    # if isnothing(datafiles)
    #     @warn "No data files will be processed!"
    # else
    #     map!(x -> joinpath(readpath, x), datafiles, datafiles)
    # end

    percent = get("percent to read", 100)

    # Qkinz
    qkinz = get("qkinz", nothing) |> abspath |> parseqkinz

    # Bad channels
    badchannels = Set{Int}(i for i in get("channels to ignore", Set{Int}()))

    atleast1d(x::AbstractArray) = x
    atleast1d(x) = [x]
    for key in keys(parameters)
        if key ∉ readkeys
            @warn "The parameter $key was supplied but not used."
        end
        if occursin("coefficient", key) || occursin("file", key)
            paths = get(key, "") |> atleast1d
            for path in paths
              if !isfile(path)
                  @warn "Could not open $key"
              end
            end
        end
    end

    Parameters(path, readpath, savepath, isotope,
               gates, calibrator,
               percent, badchannels, qkinz)
end

function make(::Type{T}, path::AbstractString, var)::T where T<:AbstractCalibrator
    if isfile(path)
        read(path, T, var=var)
    else
        T(var=var)
    end
end

function make(::Type{T}, paths::AbstractArray{<:AbstractString}, var)::T where T<:AbstractCalibrator
    base = make(T, paths[1], var)
    if length(paths) > 1
      for path in paths[2:end]
          calib = make(T, path, var)
          combine!(base, calib)
      end
    end
    base
end

### PARSE GATES ###
getgates(::Nothing) = Dict(:ex => Gate(), :t => Gate())
function getgates(unparsedgates::AbstractDict)
    # Gates on the excitation energy
    ukeys = lowercase.(keys(unparsedgates))
    ekey =  ukeys ∩ ["ex", "excitation", "ex_int"]
    promptkey = ukeys ∩ ["prompt"]
    bgkey = ukeys ∩ ["bg", "background"]

    gates = Dict{Symbol,Gate}(:ex => Gate(), :background => Gate(),
                              :promt => Gate())

    if length(ekey) == 1
        for gate in unparsedgates[ekey[1]]
            gates[:ex] = gates[:ex] ∩ parsegate(gate)
        end
    end

    if length(promptkey) == 1
        gates[:prompt] = unparsedgates[promptkey[1]] |> parsegate
    end

    if length(bgkey) == 1
        gates[:background] = unparsedgates[bgkey[1]] |> parsegate
    end

    if length(ukeys) > 0 && length(ekey ∪ promptkey ∪ bgkey) == 0
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

## Parse Qkinz
function parseqkinz(path::AbstractString)::Vector{Poly}
    dir = dirname(path)
    pattern = Regex(basename(path))
    polys = Tuple{Int, Poly{Float64}}[]
    j(x) = joinpath(dir, x)
    for fname in readdir(dir)
        m = match(pattern, fname)
        if m ≢ nothing
            f = parse(Int, m[:f])
            push!(polys, (f, JSort.BetheBloch.Banana(j(fname)).exfromede))
        end
    end
    (f, poly) = polys |> sort |> x -> zip(x...)
    @assert length(f) == 8 "Missing Qkinz strips. Got $f"
    [p for p in poly]
end
