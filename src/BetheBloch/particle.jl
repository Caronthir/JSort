include("material.jl")

struct Particle
    A::Int  # Mass Number
    Z::Int  # Element number
    material::Material  # Target material
end


function stoppingpower(p::Particle, E)
    Erel = E + p.mev
end
