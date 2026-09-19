@testset "StatisticalDistributions.LogSkewNormalFunctions" begin
    mod = ExtraStats.StatisticalDistributions.LogSkewNormalFunctions

    @test isfinite(mod.logskewnormal(1.5, 0.0, 1.0, 0.5))
    @test isfinite(mod.LogSkewNormalCDF(1.5, 0.0, 1.0, 0.5))
end
