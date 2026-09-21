export lognormal, LogNormalCDF

"""
    lognormal(x, μ, σ)

Calculate the value of a Log-Normal function with log-centre `μ` and log-standard deviation `σ` at `x`.
"""
function lognormal(x, μ, σ)
    return 1/(x*σ*sqrt(2*π)) * exp(-(log(x) - μ)^2 /(2*σ^2))
end

"""
    lognormal(x, p)

Calculate the value of a Log-Normal function with parameter vector `p` at `x`.
"""
function lognormal(x, p)
    return 1 ./(x.*p[2].*sqrt(2*π)) .* exp.(-(log.(x) .- p[1]).^2 ./(2 .*p[2].^2))
end


"""
    LogNormalCDF(x, μ, σ)

Calculate the CDF of a Log-Normal pdf with centre `μ` and standard deviation `σ` at `x`.
"""
LogNormalCDF(x, μ, σ) = ∫Normaldx((log(x) - μ) / σ)