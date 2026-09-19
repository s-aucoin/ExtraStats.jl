using ExtraStats
using Test

@testset "ExtraStats" begin
    include("MiscellaneousFunctions.jl")
    include("LeastSquaresFitting/LeastSquaresFitting.jl")
    include("StatisticalDistributions/StatisticalDistributions.jl")
end