@testset "StatisticalDistributions.LogNormalFunctions" begin
    mod = ExtraStats.StatisticalDistributions.LogNormalFunctions

    @test mod.lognormal(1.0, 0.0, 1.0) ≈ 0.3989422804014327
    @test mod.LogNormalCDF(1.0, 0.0, 1.0) ≈ 0.5
end
