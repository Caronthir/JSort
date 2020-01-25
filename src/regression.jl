using Statistics
using LinearAlgebra
using Polynomials: Poly
import PyPlot; const plt = PyPlot

function leastsquares(x, y, Ω; order::Integer = 1)
    X = ones((length(x), order+1))
    for n in 1:order
        X[:, n+1] = x.^n
    end

    if ndims(Ω) == 1
        Ω = diagm(0 => Ω)
    end

    w = inv(Ω)  # Weights defined as 1/σ²
    Xᵀ = transpose(X)
    kω = inv(Xᵀ*w*X)*Xᵀ*w
    coeffs = kω*y
    Σ = kω*Ω*kω'
    coeffs, Σ
end

function leastsquares(x, y; order::Integer = 1)
    X = ones((length(x), order+1))
    for n in 1:order
        X[:, n+1] = x.^n
    end

    Xᵀ = transpose(X)
    inv(Xᵀ*X)*Xᵀ*y
end

function plotresiduals(x, y, coefficients; ax=nothing)
    if ax ≡ nothing
        fig, ax = plt.subplots(ncols=2)
    end
    p = Poly(coefficients)
    ŷ = p.(y)
    xrange = range(minimum(x), maximum(x), length=100)
    ŷrange = p.(xrange)

    ê = ŷ - y
    ax[1].scatter(x, y)
    ax[1].scatter(x, ŷ)
    ax[1].plot(xrange, ŷrange, "--")
    ax[1].set_title("Fit")
    ax[2].scatter(x, ê)
    ax[2].set_title("Residuals")
    return ax
end
