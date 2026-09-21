using SpecialFunctions # For extra math functions

export gaussian, normal, ∫Normaldx, NormalCDF


"""
    gaussian(x, A, μ, σ)

Calculate the value of a Gaussian function with amplitude `A`, centre `μ`, and standard deviation `σ` at `x`.
"""
function gaussian(x, A, μ, σ)
    return A*exp(-(x - μ)^2/(2*σ^2))
end

"""
    gaussian(x, p)

Calculate the value of a Gaussian function with parameter vector `p` at `x`.
"""
function gaussian(x, p)
    return p[1].*exp.(-(x .- p[2]).^2 ./(2 .*p[3].^2))
end


"""
    normal(x, μ, σ)

Calculate the value of a Normal function with centre `μ` and standard deviation `σ` at `x`.
"""
function normal(x, μ, σ)
    return 1/(sqrt(2*π*σ^2)) * exp(-(x - μ)^2/(2*σ^2)) # make the pdf normalized
end

"""
    normal(x, p)

Calculate the value of a Normal function with parameter vector `p` at `x`.
"""
function normal(x, p)
    return 1/(sqrt.(2*π .* p[2].^2)) .* exp.(-(x .- p[1]).^2 ./(2 .*p[2].^2)) # make the pdf normalized
end


"""
    ∫Normaldx(x)

Calculate the value of the integral of a Normal function at `x`.
"""
∫Normaldx(x) = 1/2*(1 + erf(x / (sqrt(2))))


"""
    NormalCDF(x, μ, σ)

Calculate the CDF of a Normal pdf with centre `μ` and standard deviation `σ` at `x`.
"""
NormalCDF(x, μ, σ) = ∫Normaldx((x - μ) / σ)