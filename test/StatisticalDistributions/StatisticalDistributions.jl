@testset "StatisticalDistributions" begin
    pdf(x) = exp(-(x^2) / 2) / sqrt(2π)
    g(x) = x^2
    ginv(x) = sqrt(x)
    @test isfinite(transform_pdf(pdf, ginv)(1.0))

    @test isfinite(find_pdf_mode(pdf, 0.0))

    samples = collect(range(-3.0, 3.0, length=200))
    @test length(find_HDI(pdf, samples; target_mass=0.68)) == 2
end

include("NormalFunctions.jl")
include("SkewNormalFunctions.jl")
include("LogNormalFunctions.jl")
include("LogSkewNormalFunctions.jl")