export logskewnormal, LogSkewNormalCDF

"""
    logskewnormal(x, λ, ω, ξ)

Calculate the value of a Log-Skew-Normal function with location `λ`, scale `ω`, and shape `ξ` at `x`.
"""
function logskewnormal(x, λ, ω, ξ)
    return 2*lognormal(x, λ, ω).*∫Normaldx(ξ*(log(x) - λ)/ω)
end

"""
    logskewnormal(x, p)

Calculate the value of a Log-Skew-Normal function with parameter vector `p` at `x`.
"""
function logskewnormal(x, p)
    return 2*lognormal(x, p).*∫Normaldx(p[4]*(log.(x) - p[2])/p[3])
end


"""
    LogSkewNormalCDF(x, λ, ω, ξ)

Calculate the CDF of a Log-Skew-Normal pdf with location `λ`, scale `ω`, and shape `ξ` at `x`.
"""
LogSkewNormalCDF(x, λ, ω, ξ) = LogNormalCDF(x, λ, ω) - 2*OwensTFunction((log(x) - λ)/ω, ξ)