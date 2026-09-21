module ExtraStats

    include("MiscellaneousFunctions.jl")
    using .MiscellaneousFunctions
    export MiscellaneousFunctions

    include("LeastSquaresFitting/LeastSquaresFitting.jl")
    using .LeastSquaresFitting
    export LeastSquaresFitting

    include("StatisticalDistributions/StatisticalDistributions.jl")
    using .StatisticalDistributions
    export StatisticalDistributions

end
