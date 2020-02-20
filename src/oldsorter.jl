using .JSort
import CSV
import DelimitedFiles: writedlm
using Polynomials

const E_MIN = 5  # E Threshold

function sortevent!(matrices, event::Event, parameters::Parameters)
    #=
    Algorithm:
    For each ΔE, look for an E in the associated E-detector. If no E is found,
    or multiple correlations are made, ignore this count.
    =#
    minievent = correlate(event)

    # increment!(matrices[:emul], length(minievents))
    # increment!(matrices[:Δemul], length(minievents))

    # If there were more than 1 correlated event 
    # in the time window, it is impossible to tell which
    # LaBr corresponds to which SiRi detection
    isnothing(minievent) && return

    # Calibrate the events
    e = calibrate(minievent, :e, parameters)
    Δe = calibrate(minievent, :Δe, parameters)

    # Reject event if under threshold
    e ≤ E_MIN || Δe ≤ E_MIN && return

    # Look for events within the banana defined by Qkinz
    thick::Float64 = parameters.range(e+Δe) - parameters.range(e)
    thickdev::Float64 = parameters.thickrange[2] + e*parameters.thickrange[3]

    # Increment the matrices
    b, f = minievent.back, minievent.front
    increment!(matrices[:front], Δe, minievent.frontchannel)
    increment!(matrices[:back], e, minievent.backchannel)
    increment!(matrices[:eΔe], e, Δe)
    # increment!(matrices[Symbol("eΔeb",b, "f", f)], e, Δe)
    increment!(matrices[:thick], thick)

    # Check if the event comes within the banana
    within_banana = abs(thick-parameters.thickrange[1]) < thickdev
    !within_banana && return

    # Increment even more matrices
    increment!(matrices[:eΔethick], e, Δe)
    increment!(matrices[:heΔe], e+Δe)
    # increment!(matrices[Symbol("heΔerf", f)], e+Δe)
    # increment!(matrices[Symbol("heΔeb", b, "f", f)], e+Δe)

    ex = calibrate(e+Δe, minievent, :ex, parameters)
    increment!(matrices[:hex], ex)
    # increment!(matrices[Symbol("hexrf", f)], ex)

    # TODO Find a way to support time evolution
    # Time evolution

    # LaBr
    increment!(matrices[:labrmul], event.labrnum)
    for (channel, labr) in eachlabr(event)
        # Check that both adc and tdc have been set
        !islabrgood(labr) && return
        labre = calibrate(labr.adc, :labre, parameters)
        # TODO why this magic /8?
        labrtime = calibrate(labr.tdc/8, :labrtime, parameters)
        labrtime_c = tLaBr(labrtime, labre, e)

        increment!(matrices[:labre], labre, channel)
        increment!(matrcies[:labrtimec], channel)

    end
    # for (channel, labr) in event.labr
        # labr.adc ≤ 0 && continue
        # e = 0
    # end
    return
    # e = ecalibrator(channel_b, channel_f)(si_e[channel_b])
    # Do thick stuff here
end


function oldmakeeΔe(events::Array{TrueEvent, 1}, parameters::Parameters)
    meΔe = OMatrix{Array{Int64, 2}}((NUMBINS_E, NUMBINS_E), max=(MAX_E, MAX_ΔE))
    meΔebf = [OMatrix{Array{Int64, 2}}((NUMBINS_E, NUMBINS_E), max=(MAX_E, MAX_ΔE))
              for f in 1:8, b in 1:8]
    meΔebf_raw = [[] for f in 1:8, b in 1:8]

    for (i, event) in enumerate(events)
        e = calibrate(event.e, 0.0, 5.0)
        Δe = calibrate(event.Δe, 0.0, 2.5)
        push!(meΔebf_raw[event.b, event.f], (e, Δe))
        increment!(meΔe, e, Δe)
        increment!(meΔebf[event.b, event.f], e, Δe)
    end
    #save(meΔe, "mede", parameters.savepath)
    for b in 1:8, f in 1:8
        #save(meΔebf[b, f], "medeb$(b)f$f", parameters.savepath)
        path = joinpath(getpath(parameters, :savepath), "/medeb$(b)f$f.csv")
        CSV.write(path, meΔebf_raw[b, f])
    end
    meΔe, meΔebf
end

function makeeΔe(events::AbstractVector{TrueEvent}, parameters::Parameters; docalibrate=true)
    eΔe = [Tuple{Float32, Float32}[] for f in 1:8, b in 1:8]
    ecalibrators  = parameters.calibrator[:e]
    Δecalibrators = parameters.calibrator[:Δe]
    e = Δe = 0.0
    for (i, event) in enumerate(events)
        e  = ecalibrators[event.f, event.b](event.e)
        Δe = Δecalibrators[event.f, event.b](event.Δe)
        push!(eΔe[event.f, event.b], (e, Δe))
    end
    for f in 1:8
        total = Tuple{Float32, Float32}[]
        for b in 1:8
            path = joinpath(getpath(parameters, :savepath), "edeb$(b)f$(f).csv")
            # Stupid writedlm syntax
            data = []
            push!(data, ("e", "de"))
            append!(data, eΔe[f, b])
            open(path, "w") do io
                writedlm(io, data, ',')
            end
            append!(total, eΔe[f, b])
        end
        path = joinpath(getpath(parameters, :savepath), "edef$(f).csv")
        data = []
        push!(data, ("e", "de"))
        append!(data, total)
        open(path, "w") do io
            writedlm(io, data, ',')
        end
    end
end

function makeeΔebin(events::AbstractVector{TrueEvent}, parameters::Parameters; T=Float32)
    eΔe = [Tuple{T, T}[] for f in 1:8, b in 1:8]
    ecalibrators  = parameters.calibrator[:e]
    Δecalibrators = parameters.calibrator[:Δe]
    e = Δe = 0.0
    for (i, event) in enumerate(events)
        e  = ecalibrators[event.f, event.b](event.e)
        Δe = Δecalibrators[event.f, event.b](event.Δe)
        push!(eΔe[event.f, event.b], (e, Δe))
    end
    for f in 1:8
        path = joinpath(getpath(parameters, :savepath), "edef$(f).bin")
        iototal = open(path, "w")
        # Write the total length
        totallength = length.(eΔe[f, :]) |> sum |> T
        write(iototal, totallength)
        for b in 1:8
            path = joinpath(getpath(parameters, :savepath), "edeb$(b)f$(f).bin")
            data = eΔe[f, b]
            # Stupid writedlm syntax
            open(path, "w") do io
                write(io, length(data) |> T)
                write(io, data)
            end
            # Append to the summed data
            write(iototal, data)
        end
        close(iototal)
    end
end

function makelabr(events::AbstractVector{TrueEvent}, parameters::Parameters; docalibrate=true)
    # TODO: OMatrices are *almost* useless
    # Use long arrays instead
    #gamma_e = OMatrix{Array{Int64, 2}}((2000, 32), max=(14000, 32))
    #gamma_t = OMatrix{Array{Int64, 2}}((500, 32))
    channels = 1:33
    Γe = Dict(i => Float32[] for i in channels)
    Γt = Dict(i => Float32[] for i in channels)

    ecalibrators = parameters.calibrator[:γ]
    tcalibrators = parameters.calibrator[:t]
    @debug "Making LaBr gamma events"
    e = 0.0
    t = 0.0
    for event in events
        for γ in event.labr
            e = ecalibrators[γ.channel+1](γ.adc)
            t = tcalibrators[γ.channel+1](γ.tdc/8)
            #increment!(gamma_e, e, gamma.channel)
            push!(Γe[γ.channel+1], e)
            push!(Γt[γ.channel+1], t)
            #increment!(gamma_t, t, gamma.channel)
        end
    end
    #save(gamma_e, "gammae", getpath(parameters, :savepath))
    #save(gamma_t, "gammat", getpath(parameters, :savepath))
    return Γe, Γt
end

