using Integrals        # For computing integrals
#import ExtraStats.StatisticalDistributions.NormalFunctions: gaussian, normal, ∫Normaldx

export skewgaussian, skewnormal,
        SkewNormalCDF, OwensTFunction, SkewNormalCDF, 
        skew_δ, skew_mean, SkewNormalμ2λ,
        skew_σ, SkewNormalσ2ω,
        skew_skew, SkewNormalSkewness2Posδ, SkewNormalSkewness2Posξ,
        skew_mode

"""
    skewgaussian(x, A, λ, ω, ξ)

Calculate the value of a Skew-Gaussian function with amplitude `A`, location `λ`, scale `ω`, and shape `ξ` at `x`.
"""
function skewgaussian(x, A, λ, ω, ξ)
    return 2*gaussian(x, A, λ, ω)*∫Normaldx(ξ*(x - λ)/ω)
end

"""
    skewgaussian(x, p)

Calculate the value of a Skew-Gaussian function with parameter vector `p` at `x`.
"""
function skewgaussian(x, p)
    return 2*gaussian(x, p).*∫Normaldx.(p[4]*(x - p[2])/p[3])
end


"""
    skewnormal(x, λ, ω, ξ)

Calculate the value of a Skew-Normal function with location `λ`, scale `ω`, and shape `ξ` at `x`.
"""
    function skewnormal(x, λ, ω, ξ)
        return 2*normal(x, λ, ω)*∫Normaldx(ξ*(x - λ)/ω)
    end

"""
    skewnormal(x, p)

Calculate the value of a Skew-Normal function with parameter vector `p` at `x`.
"""
    function skewnormal(x, p)
        return 2*normal(x, p).*∫Normaldx(p[3].*(x .- p[1])./p[2])
    end


"""
    SkewNormalCDF(x, λ, ω, ξ)

Calculate the CDF of a Skew-Normal pdf with location `λ`, scale `ω`, and shape `ξ` at `x`.
"""
    SkewNormalCDF(x, λ, ω, ξ) = ∫Normaldx((x - λ) / ω) - 2*OwensTFunction((x - λ)/ω, ξ)


"""
    OwensTFunction(h, a)

Calculate the value of Owen's T function with parameters `h` and `a`.
"""
function OwensTFunction(h, a)
    integrand(x, p) = 1/(1 + x^2) * exp(-1/2 * p^2 * (1 + x^2))
    bounds = (0, a)
    int_prob = IntegralProblem(integrand, bounds, h)
    return 1/(2π) * solve(int_prob, QuadGKJL()).u
end


###################################
# Skew-Normal Statistics functions #

"""
    skew_δ(ξ)

Calculate the convienence function `δ` from shape `ξ`.
"""
function skew_δ(ξ)
    return ξ / sqrt(1 + ξ^2)
end


"""
    skew_mean(λ, ω, ξ)

Calculate the mean of a Skew-Normal function with location `λ`, scale `ω`, and shape `ξ`.
"""
function skew_mean(λ, ω, ξ)
    return λ + ω*skew_δ(ξ)*sqrt(2/π)
end


"""
    SkewNormalμ2λ(μ, σ, k)

Calculate the location `λ` of a Skew-Normal function with mean `μ`, standard deviation `σ`, and skewness `k`.
"""
SkewNormalμ2λ(μ, σ, k) = μ - sqrt(2/π) * SkewNormalσ2ω(σ, k)


"""
    skew_σ(ω, ξ)

Calculate the standard deviation of a Skew-Normal function with scale `ω` and shape `ξ`.
"""
function skew_σ(ω, ξ)
    return ω * sqrt(1 - (2 * skew_δ(ξ)^2 / π))
end


"""
    SkewNormalσ2ω(μ, σ, k)

Calculate the scale parameter `ω` of a Skew-Normal function with standard deviation `σ` and skewness `k`.
"""
SkewNormalσ2ω(σ, k) = σ * (1 - 2/π * SkewNormalSkewness2Posδ(k)^2)^(-1/2)


"""
    skew_skew(ξ)

Calculate the skewness of a Skew-Normal function with shape `ξ`.
"""
function skew_skew(ξ)
    return (4 - π)/2 * (skew_δ(ξ)*sqrt(2/π))^3 / (1 - 2*skew_δ(ξ)/π)^(3/2)
end


"""
    SkewNormalSkewness2Posδ(k)

Calculate the convienence function `δ` of a Skew-Normal function with skewness `k`.
"""
function SkewNormalSkewness2Posδ(k)
    # I think this equation is only valid for positive shape parameters/skewness values
    radical(k) = (4 - π)^(1/3)*k^(4/3) + π*(2)^(1/3)*(4 - π)*k^(2/3) # equation for the term in the square root
    return (sqrt(radical(k)) - (4 - π)^(1/6)*k^(2/3)) / ((2)^(1/3) * (4 - π)^(5/6)) # positive root of the solution for δ
end


"""
    SkewNormalSkewness2Posξ(k)

Calculate the shape parameter `ξ` of a Skew-Normal function with skewness `k`.
"""
function SkewNormalSkewness2Posξ(k)
    return SkewNormalSkewness2Posδ(k) / sqrt(1 - SkewNormalSkewness2Posδ(k)^2) # convert δ to ξ
end


## Mode convience function ###
m₀(ξ) = sqrt(2/π)*skew_δ(ξ) - (1 - π/4)*(sqrt(2/π)*skew_δ(ξ))^3/(1 - 2/π*skew_δ(ξ)^2) - sign(ξ)/2*exp(-2*π/abs(ξ))

"""
    skew_mode(λ, ω, ξ)

Calculate the mode of a Skew-Normal function with location `λ`, scale `ω`, and shape `ξ`.
"""
function skew_mode(λ, ω, ξ)
    return λ + ω*m₀(ξ)
end