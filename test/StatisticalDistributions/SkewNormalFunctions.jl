@testset "StatisticalDistributions.SkewNormalFunctions" begin
    @test isfinite(skewgaussian(0.0, 1.0, 0.0, 1.0, 0.5))
    @test isfinite(skewnormal(0.0, 0.0, 1.0, 0.5))
    @test isfinite(OwensTFunction(0.0, 0.5))
    @test isfinite(skew_δ(0.5))
    @test isfinite(skew_mean(0.0, 1.0, 0.5))
    @test isfinite(skew_σ(1.0, 0.5))
    @test isfinite(skew_skew(0.5))
    @test isfinite(skew_mode(0.0, 1.0, 0.5))
end
