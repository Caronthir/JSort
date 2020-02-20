using Polynomials: Poly

struct Banana{T}
    ex   ::Vector{T}
    Δe   ::Vector{T}
    dΔe  ::Vector{T}
    e    ::Vector{T}
    de   ::Vector{T}
    etot ::Vector{T}
    exfromede::Poly
end

function Banana(filename::Union{String, IO})
    rows = Vector{Float32}[]
    poly = nothing
    for line in readlines(filename)
        tokens = split(line, r",?\s+")
        if length(tokens) != 6
            if length(tokens) > 1 && startswith(tokens[1], "chi")
                a0 = tokens[5][1:end-3]
                a1 = tokens[8]
                a2 = tokens[11][1:end-8]
                (a0, a1, a2) = parse.(Float32, [a0, a1, a2])
                poly = Poly{Float32}([a0*1e3, a1, a2/1e3])
            end
            continue
        end
        try
          push!(rows, parse.(Float32, tokens))
        catch ArgumentError
            continue
        end
    end
    rows = hcat(rows...)
    Banana{Float32}([rows[i, :] for i in 1:6]..., poly)
end

function smash(banana::Banana)

end
