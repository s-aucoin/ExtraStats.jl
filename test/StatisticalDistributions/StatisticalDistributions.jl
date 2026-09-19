@testset "StatisticalDistributions" begin
    mod = ExtraStats.StatisticalDistributions

    pdf(x) = exp(-(x^2) / 2) / sqrt(2π)
    g(x) = x^2
    ginv(x) = sqrt(x)
    @test isfinite(mod.transform_pdf(pdf, ginv)(1.0))

    @test isfinite(mod.find_pdf_mode(pdf, 0.0))

    samples = collect(range(-3.0, 3.0, length=200))
    @test length(mod.find_HDI(pdf, samples; target_mass=0.68)) == 2
end

include("NormalFunctions.jl")
include("SkewNormalFunctions.jl")
include("LogNormalFunctions.jl")
include("LogSkewNormalFunctions.jl")