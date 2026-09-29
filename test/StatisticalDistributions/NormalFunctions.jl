@testset "StatisticalDistributions.NormalFunctions" begin
    @test gaussian(0.0, 1.0, 0.0, 1.0) ≈ 1.0
    @test normal(0.0, 0.0, 1.0) ≈ 0.3989422804014327
    @test ∫Normaldx(0.0) ≈ 0.5
    @test NormalCDF(0.0, 0.0, 1.0) ≈ 0.5
end
