@testset "StatisticalDistributions.LogSkewNormalFunctions" begin
    @test isfinite(logskewnormal(1.5, 0.0, 1.0, 0.5))
    @test isfinite(LogSkewNormalCDF(1.5, 0.0, 1.0, 0.5))
end
