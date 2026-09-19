@testset "LeastSquaresFitting.LinearFitUncertainty" begin
    mod = ExtraStats.LeastSquaresFitting.LinearFitUncertainty

    xi = [1.0, 2.0, 3.0, 4.0]
    yi = [3.0, 5.0, 7.0, 9.0]
    yhat = 2.0 .* xi .+ 1.0

    ci = mod.linear_uncertainty(xi, yi, yhat, 0.95)
    @test isfinite(ci.beta_e)
    @test isfinite(ci.alpha_e)

    @test isfinite(mod.SE_beta(xi, yi, yhat))
    @test isfinite(mod.SE_alpha(xi, yi, yhat))

    xvals = collect(0.0:0.5:5.0)
    lines = mod.uncertainty_lines(2.0, 1.0, xvals, xi, yi, yhat, 0.95)
    @test length(lines.upper_ln) == length(xvals)
    @test length(lines.lower_ln) == length(xvals)
end
