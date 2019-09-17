using .JSort
import CSV

const E_MIN = 5  # E Threshold

function correlate(event::Event)::Union{Nothing, MiniEvent}
    foundevent = false
    let minievent
        for Δei in 1:event.Δenum
            Δe = event.Δe[Δei]
            for ei in 1:event.enum
                e = event.e[ei]
                back = e.channel
                if Δe.associated_channel ≠ back
                    continue
                end
                front = Δe.channel % 8
                # We can only handle 1 correlated event
                if foundevent
                    return nothing
                end 
                foundevent = true
                minievent = MiniEvent(front+1, back+1, Δe.channel, e.channel, Δe.data, e.data)
            end
        end
        if foundevent
            return minievent
        else
            return nothing
        end
    end
end


function sortevent!(matrices, event::Event, parameters::Parameters)
    #=
    Algorithm:
    For each ΔE, look for an E in the associated E-detector. If no E is found,
    or multiple corelations are made, ignore this count.
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

function correlateevent(event::Event)
    minievent = correlate(event)
    isnothing(minievent) && return
    labrs = LaBrEvent(minievent)
    isgood = false
    for (channel, labr) in eachlabr(event)
        if islabrgood(labr)
            isgood = true
            push!(labrs.labr, ADCTDC(labr.channel, labr.adc, labr.tdc))
        end
    end
    !isgood && return
    return labrs
end

function makeeΔe(events::Array{TrueEvent, 1}, parameters::Parameters)
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
        CSV.write(parameters.savepath * "/medeb$(b)f$f.csv", meΔebf_raw[b, f])
    end
    meΔe, meΔebf
end

function makelabr(events::Array{TrueEvent, 1}, parameters::Parameters)
    gamma_e = OMatrix{Array{Int64, 2}}((2000, 32), max=(14000, 32))
    gamma_t = OMatrix{Array{Int64, 2}}((500, 32))
    for event in events
        for gamma in event.labr
            e = calibrate(gamma.adc, 0, 5)
            t = calibrate(gamma.tdc/8, 0, 1)
            increment!(gamma_e, e, gamma.channel)
            increment!(gamma_t, t, gamma.channel)
        end
    end
    save(gamma_e, "gammae", parameters.savepath)
    save(gamma_t, "gammat", parameters.savepath)
    return gamma_e, gamma_t
end
