using Documenter
using ExtraStats

makedocs(
    sitename = "ExtraStats",
    format = Documenter.HTML(),
    modules = [ExtraStats],
    remotes = nothing,
    pages = ["Library" => "library.md", 
            "Miscellaneous" => "misc.md", 
            "Least Squares Fitting" => "LeastSquaresFitting/least_squares_fitting.md", 
            "Statistical Distributions" => "StatisticalDistributions/statistical_distributions.md"]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
