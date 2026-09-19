module FitEvaluation

using Statistics
using StatsBase        # For more stats
using Optim            # For optimizing functions
using Distributions    # For statistical distributions
import ExtraStats.MiscellaneousFunctions: nanmean, bincounts

export CalculateRMSE, r², χ²Test, KuiperTest


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
    KuiperTest(x, Fₕ)

Perform Kuiper's goodness-of-fit test of the null hypothesis that the data `x` come from CDF `Fₕ`.
"""
function KuiperTest(x, Fₕ)

    F = StatsBase.ecdf(x) # calculate the empirical CDF of the data

    PosSide(x) = F(x) - Fₕ(x) # for the positive side of the Kuiper statistic
    NegSide(x) = Fₕ(x) - F(x) # for the negative side of the Kuiper statistic

    # Find the maximum of each side using optimization #
    posopt = optimize(x -> -PosSide(x), minimum(x), maximum(x))
    PosSideMax = -Optim.minimum(posopt)

    negopt = optimize(x -> -NegSide(x), minimum(x), maximum(x))
    NegSideMax = -Optim.minimum(negopt)


    return PosSideMax + NegSideMax # calculate the Kuiper statistic
end

end