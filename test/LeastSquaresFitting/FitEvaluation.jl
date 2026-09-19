using Distributions

@testset "LeastSquaresFitting.FitEvaluation" begin
    mod = ExtraStats.LeastSquaresFitting.FitEvaluation

    @test mod.CalculateRMSE([3.0, 4.0]) ≈ sqrt(12.5)

    data = [1.0, 2.0, 3.0, 4.0]
    resid = [0.0, 0.0, 0.0, 0.0]
    @test mod.r²(data, resid) ≈ 1.0

    x = randn(100)
    dist_cdf(z, p) = map(Base.Fix1(cdf, Normal(p[1], p[2])), z)
    res = mod.χ²Test(x, dist_cdf, [0.0, 1.0]; nbins=20, siglevel=0.95)
    @test isfinite(res.χ²_stat)
    @test isfinite(res.χ²_crit)

    Fh = z -> map(Base.Fix1(cdf, Normal(0.0, 1.0)), z)
    @test mod.KuiperTest(x, Fh) >= 0.0
end
