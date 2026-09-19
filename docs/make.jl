using Documenter
using ExtraStats

makedocs(
    sitename = "ExtraStats",
    format = Documenter.HTML(),
    modules = [ExtraStats]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
