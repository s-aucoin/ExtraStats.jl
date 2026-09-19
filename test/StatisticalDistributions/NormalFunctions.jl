@testset "StatisticalDistributions.NormalFunctions" begin
    mod = ExtraStats.StatisticalDistributions.NormalFunctions

    @test mod.gaussian(0.0, 1.0, 0.0, 1.0) ≈ 1.0
    @test mod.normal(0.0, 0.0, 1.0) ≈ 0.3989422804014327
    @test mod.∫Normaldx(0.0) ≈ 0.5
    @test mod.NormalCDF(0.0, 0.0, 1.0) ≈ 0.5
end
