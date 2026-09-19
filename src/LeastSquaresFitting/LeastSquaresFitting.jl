module LeastSquaresFitting

using StatsBase        # For more stats
using HypothesisTests  # For even more stats
using LinearAlgebra    # For linear algebra
using LsqFit           # For least squares fitting
using CairoMakie       # For plotting

export fittest, OLS, GLS, sqcovar

include("FitEvaluation.jl")
using .FitEvaluation

include("LinearFitUncertainty.jl")
using .LinearFitUncertainty


"""
    fittest(x, y; fittype = :OLS, fitmodel = nothing, p0 = [0.1, 0.1], figshow = true)

Perform a least squares fit with dependant variable `y` and independant `x`, test the significance of the trend,
and create a plot to evaluate the residuals.

Default `fitmodel` is linear, and `fittype` ordinary least squares fit, `:OLS`, but optionally perform generalized least squares with `fittype = :GLS`.

# Example
```
julia> using ExtraStats, Random

julia> Random.seed!(1234)
TaskLocalRNG()

julia> x = 0.1:0.1:100
0.1:0.1:100.0

julia> y = randn(1000)
1000-element Vector{Float64}:
  0.9706563288552144
  0.871497852880908
  ⋮
  0.6530055215854716
 -1.2911234902314417

julia> fmodel(t, p) = p[1] * t .+ p[2]
fmodel (generic function with 1 method)

julia> fittest(x, y; fittype = :OLS, fitmodel = fmodel, p0 = [0.1, 0.1], figshow = true)
(LsqFit.LsqFitResult{Vector{Float64}, [...]}
------   --------------   --------------
, Float64[]), Test for nonzero correlation
----------------------------
Population details:
    parameter of interest:   Correlation
    value under h_0:         0.0
    point estimate:          0.0611079
    95% confidence interval: (-0.0008886, 0.1226)

Test summary:
    outcome with 95% confidence: fail to reject h_0
    two-sided p-value:           0.0534

Details:
    number of observations:          1000
    number of conditional variables: 0
    t-statistic:                     1.93408
    degrees of freedom:              998
, Scene (1200px, 750px):
  0 Plots
  15 Child Scenes:
    ├ Scene (1200px, 750px)
    [...]
```
"""
function fittest(x, y; fittype = :OLS, fitmodel = nothing, p0 = [0.1, 0.1], figshow = true)

    fitfunc = getfield(LeastSquaresFitting, fittype)

    LSfit = fitfunc(x, y; fitmodel = fitmodel, p0 = p0)


    ### Test the correlation significance ###
    ttest = CorrelationTest(x, y)


    if figshow
        ### Make the Plot ###
        ## Define the kwargs ##
        xlabgap = 30
        ylabgap = 50
        axis_kwargs = ()
        plot_kwargs = ()
        label_kwargs = (fontsize=24, font=:bold)

        ## Initialize the figure ##
        fig = Figure(size=(1200, 750))

        ## Create the axes ##
        ax1 = Axis(fig[1,1]; axis_kwargs...)
        ax2 = Axis(fig[2,1]; axis_kwargs...)
        ax3 = Axis(fig[1,2]; axis_kwargs...)
        ax4 = Axis(fig[2,2]; axis_kwargs...)

        ## Make the plots ##
        # Input scatter and fit line
        ax1_xlab = "x"
        ax1_ylab = "y"
        # Use the matrix K to transform the data back
        scatter!(ax1, x, y)
        lines!(ax1, x, fitmodel(x, LSfit.param), color=:orange, label="Fit")
        axislegend(ax1, position=:rt, framevisible=false, labelsize=16)
        Label(fig[1, 1, Bottom()], ax1_xlab; padding=(0,0,0,xlabgap), label_kwargs...)
        Label(fig[1, 1, Left()], ax1_ylab; padding=(0,ylabgap,0,0), rotation = pi/2, label_kwargs...)


        # Residuals vs input x
        ax2_xlab = "x"
        ax2_ylab = "Residuals r"
        scatter!(ax2, x, LSfit.resid)
        Label(fig[:, 1, Bottom()], ax2_xlab; padding=(0,0,0,xlabgap), label_kwargs...)
        Label(fig[2, 1, Left()], ax2_ylab; padding=(0,ylabgap,0,0), rotation = pi/2, label_kwargs...)


        # Residual scatter with itself shifted by 1
        ax3_xlab = "r(n)"
        ax3_ylab = "r(n+1)"
        scatter!(ax3, LSfit.resid[1:end-1], LSfit.resid[2:end])
        Label(fig[1, 2, Bottom()], ax3_xlab; padding=(0,0,0,xlabgap), label_kwargs...)
        Label(fig[1, 2, Left()], ax3_ylab; padding=(0,ylabgap,0,0), rotation = pi/2, label_kwargs...)


        # Autocorrelation of the residuals
        ax4_xlab = "Lag τ"
        ax4_ylab = "ρᵣᵣ"
        lows = zeros(length(autocor(LSfit.resid)))
        rangebars!(ax4, 1:length(autocor(LSfit.resid)), lows, autocor(LSfit.resid))

        # 95% confidence for gaussian noise
        conf = sqrt(2)*erfinv(0.95) / sqrt(length(LSfit.resid))
        hlines!(ax4, conf, linestyle=:dash, color=:orange, label="95% Confidence")
        hlines!(ax4, -conf, linestyle=:dash, color=:orange)
        axislegend(ax4, position=:rt, framevisible=false, labelsize=16)

        Label(fig[2, 2, Bottom()], ax4_xlab; padding=(0,0,0,xlabgap), label_kwargs...)
        Label(fig[2, 2, Left()], ax4_ylab; padding=(0,ylabgap,0,0), rotation = pi/2, label_kwargs...)


        colgap!(fig.layout, 20)
        rowgap!(fig.layout, 20)

        ## Global Title ##
        title = "Fit Test, p-value = $(pvalue(ttest))"
        Label(fig[:, :, Top()], title; padding=(0,0,20,0), label_kwargs...)

    else
        fig = nothing
    end

    return (; LSfit, ttest, fig)
end



"""
    OLS(x, y; fitmodel = nothing, p0 = [0.1, 0.1])

Perform an ordinary least squares fit with dependant variable `y` and independant `x`.

Default `fitmodel` is linear.

# Example
```
julia> using ExtraStats, Random

julia> Random.seed!(1234)
TaskLocalRNG()

julia> x = 0.1:0.1:100
0.1:0.1:100.0

julia> y = randn(1000)
1000-element Vector{Float64}:
  0.9706563288552144
  0.871497852880908
  ⋮
  0.6530055215854716
 -1.2911234902314417

julia> fmodel(t, p) = p[1] * t .+ p[2]
fmodel (generic function with 1 method)

julia> OLS(x, y; fitmodel = fmodel, p0 = [0.1, 0.1])
LsqFit.LsqFitResult{Vector{Float64}, [...]}
```
"""
function OLS(x, y; fitmodel = nothing, p0 = [0.1, 0.1])
    # By default perform a linear fit
    if fitmodel == nothing
        fmodel(t, p) = p[1] * t .+ p[2]
        fitmodel = fmodel
    end

    return curve_fit(fitmodel, x, y, p0)
end


"""
    GLS(x, y; K = nothing, fitmodel = nothing, p0 = [0.1, 0.1])

Perform a generalized least squares fit with dependant variable `y` and independant `x`.

Default `fitmodel` is linear.

Weight PrecisionMatrix `K` is estimated from the residuals of OLS by default.

# Example
```
julia> using ExtraStats, Random

julia> Random.seed!(1234)
TaskLocalRNG()

julia> x = 0.1:0.1:100
0.1:0.1:100.0

julia> y = randn(1000)
1000-element Vector{Float64}:
  0.9706563288552144
  0.871497852880908
  ⋮
  0.6530055215854716
 -1.2911234902314417

julia> fmodel(t, p) = p[1] * t .+ p[2]
fmodel (generic function with 1 method)

julia> GLS(x, y; fitmodel = fmodel, p0 = [0.1, 0.1])
LsqFit.LsqFitResult{Vector{Float64}, [...]}
```
"""
function GLS(x, y; K = nothing, fitmodel = nothing, p0 = [0.1, 0.1])
    # By default perform a linear fit
    if fitmodel == nothing
        fmodel(t, p) = p[1] * t .+ p[2]
        fitmodel = fmodel
    end

    if isnothing(K)
        OLS_fit = OLS(x, y; fitmodel = fitmodel, p0 = p0)
        K = PrecisionMatrix(inv(sqcovar(OLS_fit.resid)))
    end

    return curve_fit(fitmodel, x, y, K, p0)
end



"""
    sqcovar(X)

Calculate the Hermitian square root auto covariance matrix of input vector `X`.
"""
function sqcovar(X)

    n = length(X)
    rho = autocor(X, collect(0:n-1))

    V = Matrix{Float64}(undef, n, n)
    for (i, j) in Iterators.product(1:n, 1:n)
        V[i, j] = rho[abs(i - j) + 1]
    end

    eigval, eigvec = eigen(Hermitian(V))

    K = eigvec * Diagonal(sqrt.(max.(eigval, 0.0))) * inv(eigvec)

    return Hermitian(K)
end


end