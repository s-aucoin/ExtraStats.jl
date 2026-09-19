@testset "StatisticalDistributions.SkewNormalFunctions" begin
    mod = ExtraStats.StatisticalDistributions.SkewNormalFunctions

    @test isfinite(mod.skewgaussian(0.0, 1.0, 0.0, 1.0, 0.5))
    @test isfinite(mod.skewnormal(0.0, 0.0, 1.0, 0.5))
    @test isfinite(mod.OwensTFunction(0.0, 0.5))
    @test isfinite(mod.skew_δ(0.5))
    @test isfinite(mod.skew_mean(0.0, 1.0, 0.5))
    @test isfinite(mod.skew_σ(1.0, 0.5))
    @test isfinite(mod.skew_skew(0.5))
    @test isfinite(mod.skew_mode(0.0, 1.0, 0.5))
end
