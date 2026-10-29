using Statistics
using StatsBase        # For more stats
using Distributions    # For statistical distributions

export CalculateRMSE, r², χ²Test, KuiperTest, KuiperStat, KuiperFPP, KSTest, KSStat, KSPValue


"""
    CalculateRMSE(e)

Calculate the Root-Mean-Square-Error of the vector of residuals `e`.
"""
CalculateRMSE(e) = sqrt(mean(e.^2))


"""
    r²(data, e)

Estimate the coefficient of determination of a fit to `data` with residuals `e`.
"""
function r²(data, e)
    SSR = sum(e.^2)
    SST = sum((filter(!isnan, data) .- nanmean(data)).^2)
    return 1 - SSR/SST
end


"""
    χ²Test(x, CDF, p; nbins=500, siglevel = 0.99)

Perform a χ² goodness-of-fit test of the null hypothesis that the data `x` come from a distribution with cumulative distribution function `cdf` with parameters `p`.

Optionally specify the number of bins `nbins` to separate `x` into, and the significance level `siglevel` to test against.
"""
function χ²Test(x, CDF, p; nbins=500, siglevel = 0.99)

    binedges = range(minimum(x), maximum(x), length = nbins) # separate the data into bins

    Oi = bincounts(x, binedges; makepdf=false).counts # count the number of observations in each bin

    # Degrees of freedom #
    k = length(Oi)
    c = length(p) + 1

    χ²_crit = quantile(Chisq(k-c), 1 - siglevel) # calculate the χ² critical value

    discrete_CDF = CDF(binedges, p) # evaluate the CDF at the bin edges
    Ei = length(x)*(discrete_CDF[2:end] - discrete_CDF[1:end-1]) # turn into the number of expected observations in each bin

    χ²_stat = sum((Oi .- Ei).^2 ./ Ei) # calculate the χ² statistic

    if χ²_stat > χ²_crit
        println("The null hypothesis can be rejected at the $(siglevel*100)% significance level.")
    else
        println("The null hypothesis cannot be rejected at the $(siglevel*100)% significance level.")
    end

    return (; Oi, Ei, χ²_stat, χ²_crit)
end



"""
    KuiperTest(x::AbstractVector, Fₕ::Function)

Calculate the Kuiper statistic and corresponding False-Positive Probability (FPP) for observations `x` and CDF `Fₕ`.
"""
function KuiperTest(x::AbstractVector, Fₕ::Function)
    V = KuiperStat(x, Fₕ)
    N = length(x)
    P = KuiperFPP(V, N)

    return (; V, P)
end


"""
    KuiperStat(x::AbstractVector, Fₕ::Function)

Calculate the Kuiper statistic for the goodness-of-fit test of the null hypothesis that the data `x` come from CDF `Fₕ`.
"""
function KuiperStat(x::AbstractVector, Fₕ::Function)
    F = StatsBase.ecdf(x) # calculate the empirical CDF of the data

    Δ(x) = F(x) - Fₕ(x) # difference between the empirical CDF and the hypothesized CDF

    PosSideMax = maximum(Δ.(x))  # the maximum of the positive side
    NegSideMax = maximum(-Δ.(x)) # the maximum of the negative side

    return PosSideMax + NegSideMax # calculate the Kuiper statistic
end


"""
    KuiperFPP(V::AbstractFloat, N::Integer)

Calculate the False Positive Probability (FPP) given a Kuiper statistic `V` and sample size `N`.

Stephens (Biometrika, 1965) give a set of exact formulae for the False Positive Probability (FPP) valid for certain ranges of Ν (cases 1-3)
in addition to the asymptotic formulae for large Ν given by Kuiper (Indagationes Mathematicae (Proceedings), 1960) (case 4).

The interpretation of the FPP is that it is the probability that the Kuiper statistic D is larger than a given `V` under the null hypothesis,
i.e. the probability that the observed data come from the hypothesis PDF.
The null hypothesis is that the set of {xᵢ}, i=1,...,N is an outcome of `N` draws of the hypothesis PDF fₕ.

This function was written with reference to the astropy package in python.
"""
function KuiperFPP(V::AbstractFloat, N::Integer)
    VN = V*N

    if VN < 1 || V > 2
        throw(DomainError(V, "1/N ≤ V ≤ 2 by definition. (VN = $VN)"))
    end

    # (case 1)
    if VN <= 2
        return kuiperfpp_case1(V, N) 

    # (case 2)
    elseif VN <= 3
        return kuiperfpp_case2(V, N)

    # (case 3)
    elseif (V >= 1/2 && N % 2 == 0) || (VN >= (N-1)/2 && N % 2 == 1)
        return kuiperfpp_case3(V, N)

    # (case 4)
    else
        return kuiperfpp_case4(V, N)
    end
end

function kuiperfpp_case1(V::AbstractFloat, N::Integer) 
    return 1 - factorial(N) * (V - 1/N)^(N-1)
end

function kuiperfpp_case2(V::AbstractFloat, N::Integer)
    s = 1/2*sqrt(-N^2*V^2 + 6N*V - 7)
    x = 1/2*(N*V - 1)
    α = x + s
    β = x - s

    return 1 - factorial(N - 1) * (β^(N-1)*(1-α) - α^(N-1)*(1-β)) / (N^(N-2)*(β-α))
end

function kuiperfpp_case3(V::AbstractFloat, N::Integer)
    y(t) = V + t/N
    Tt(t) = y(t)^(t-3) * (N*y(t)^3 - y(t)^2*t*(3-2/N)/N - t*(t-1)(t-2)/N^2)
    summand(t) = binomial(N, t) * (1 - V - t/N)^(N-t-1) * Tt(t)

    ts = 1:(floor(Int, N * (1 - V)))
    return sum(summand.(ts))
end

function kuiperfpp_case4(V::AbstractFloat, N::Integer)
    z = V*sqrt(N)
    s1(m) = 2 * (4 * m^2*z^2 - 1)*exp(-2 * m^2*z^2)
    s2(m) = m^2*(4 * m^2*z^2 - 3)*exp(-2 * m^2*z^2)

    ms = 1:ceil(Int, 18.82 / z) # When m*z>18.82, exp(-2 * m^2*z^2) underflows so there is no point in going higher
    return sum(s1.(ms) .- 8*V/3 .* s2.(ms))
end


"""
    KSTest(x::AbstractVector, Fₕ::Function)

Calculate the Kolmogorov-Smirnov statistic and corresponding False-Positive Probability (FPP) for observations `x` and CDF `Fₕ`.
"""
function KSTest(x::AbstractVector, Fₕ::Function) 
    K = KSStat(x, Fₕ)
    N = length(x)
    P = KSPValue(K, N)

    return (K, P)
end

"""
    KSStat(x::AbstractVector, Fₕ::Function)

Calculate the Kolmogorov-Smirnov statistic for observations `x`, hypothesized to come from CDF `Fₕ`.
"""
function KSStat(x::AbstractVector, Fₕ::Function)
    F = StatsBase.ecdf(x) # calculate the empirical CDF of the data

    Δ(x) = F(x) - Fₕ(x) # difference between the empirical CDF and the hypothesized CDF
    Δ2(x)= Fₕ(x) - F(x) # difference between the empirical CDF and the hypothesized CDF

    δp = maximum(Δ.(x)) # maximum positive difference
    δn = maximum(Δ2.(x)) # maximum negative difference

    return maximum([δp, δn])
end

"""
    KSPValue(K::AbstractFloat, N::Integer)

Calculate the False Positive Probability (FPP) given a Kolmogorov-Smirnov statistic `K` and sample size `N`.

The interpretation of the FPP is that it is the probability that the K-S statistic is larger than a given `K` under the null hypothesis,
i.e. the probability that the observed data come from the hypothesis PDF.
The null hypothesis H₀ is that the set of {xᵢ}, i=1,...,N is an outcome of `N` draws of the hypothesis PDF fₕ.
The equivalent two-sided p-value is 2FPP, so e.g. at the 95% confidence level, reject H₀ if P < 0.025.
"""
function KSPValue(K::AbstractFloat, N::Integer)
    x = K * sqrt(N)
    summand(k) = (-1)^(k-1) * exp(-2 * k^2*x^2)

    ks = 1:ceil(Int, 18.82 / x)
    return sum(summand.(ks))
end