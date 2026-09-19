@testset "LeastSquaresFitting" begin
    mod = ExtraStats.LeastSquaresFitting

    x = collect(0.0:0.5:4.0)
    y = 2.0 .* x .+ 1.0
    y_obs = y .+ 0.05 .* randn(length(y))
    fmodel(t, p) = p[1] .* t .+ p[2]

    fit = mod.OLS(x, y; fitmodel=fmodel, p0=[1.0, 1.0])
    @test fit.param[1] ≈ 2.0 atol = 1e-3
    @test fit.param[2] ≈ 1.0 atol = 1e-3

    fit_gls = mod.GLS(x, y_obs; fitmodel=fmodel, p0=[1.0, 1.0])
    @test fit_gls.param[1] ≈ 2.0 atol = 0.1
    @test fit_gls.param[2] ≈ 1.0 atol = 0.1

    cov = mod.sqcovar(y_obs)
    @test size(cov) == (length(y_obs), length(y_obs))
    @test isfinite(cov[1, 1])
end

include("FitEvaluation.jl")
include("LinearFitUncertainty.jl")