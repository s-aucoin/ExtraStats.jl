@testset "StatisticalDistributions.LogNormalFunctions" begin
    @test lognormal(1.0, 0.0, 1.0) ≈ 0.3989422804014327
    @test LogNormalCDF(1.0, 0.0, 1.0) ≈ 0.5
end
