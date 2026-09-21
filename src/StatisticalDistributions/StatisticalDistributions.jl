using Symbolics        # For symbolic math
using Optim            # For optimizing functions
using NonlinearSolve   # For solving nonlinear equations
import ADTypes          # For modern autodiff backend selection

export find_pdf_mode, find_HDI, transform_pdf

include("NormalFunctions.jl")
include("LogNormalFunctions.jl")
include("SkewNormalFunctions.jl")
include("LogSkewNormalFunctions.jl")


"""
    transform_pdf(pdf, g, g⁻¹)

Complete a change of variables on an analytic probability density function `pdf` using the inverse `g⁻¹` of the transformation function.
"""
function transform_pdf(pdf, g⁻¹)
    @variables y
    dg⁻¹dy = build_function(Symbolics.derivative(g⁻¹(y), y), y; expression=Val{false})

    return build_function(pdf(g⁻¹(y)) * dg⁻¹dy(y), y; expression=Val{false})
end


"""
    find_pdf_mode(pdf, x0; lr=nothing, up=nothing)

Estimate the mode of a probability density function `pdf` with initial guess `x0`.
"""
function find_pdf_mode(pdf, x0; lr=nothing, up=nothing)

    f(x) = -pdf(x[1])
    ad = ADTypes.AutoForwardDiff()

    if lr == nothing && up == nothing
        opt = optimize(f, [x0], LBFGS(); autodiff = ad)

    elseif lr != nothing && up == nothing
        up = Inf
        opt = optimize(f, [lr], [up], [x0], Fminbox(LBFGS()); autodiff = ad)

    elseif lr == nothing && up != nothing
        lr = -Inf
        opt = optimize(f, [lr], [up], [x0], Fminbox(LBFGS()); autodiff = ad)

    elseif lr != nothing && up != nothing
        opt = optimize(f, [lr], [up], [x0], Fminbox(LBFGS()); autodiff = ad)
    end

    return Optim.minimizer(opt)[1]
end


"""
    find_HDI(pdf, samples; target_mass = 0.68)

Find the Highest Density Interval (Hyndman 1996) for the analytic probability density function `pdf`.
"""
function find_HDI(pdf, samples; target_mass = 0.68, lrguess=0.01*minimum(samples), upguess=maximum(samples))

    N = length(samples) # number of samples from pdf

    Y = pdf.(samples) # put the samples back into the pdf (CRazy right?).
    # The percentiles of this quantity are the values fα for which the HDR is defined

    Y_sorted = sort(Y, rev = true) # sort the transformed samples in descending order

    target_idx = floor(Int, target_mass*N) # the value of the percentile is related to the index of sorted values

    fα = Y_sorted[target_idx] # the probability density in pdf for which P(pdf >= fα) = 1 - α


    ## Now find the x values for where the pdf and fα intersect to define the HDI ##
    HDI_roots_f(x, p) = pdf.(x) .- fα

    ## get the initial guesses ##
    # idea is to approximate where the function changes sign #
    n = 100000 # number of points to sample
    xs = range(lrguess, upguess, length=n) # points to sample
    samp_HDI = HDI_roots_f.(xs, 0) # approximated function
    root_guess_idx = findall(samp_HDI[1:end-1] .* samp_HDI[2:end] .< 0) # Where this quantity is negative is where the function crosses 0

    if size(root_guess_idx, 1) < 2 # incase it can't find the lower crossing because of resolution problems
        root_guess_idx = [1, root_guess_idx[1]]
    end

    root_guess = xs[root_guess_idx] # make the guess the negative values of the above

    HDI_prob = NonlinearProblem(HDI_roots_f, root_guess) # set up the root finding problem
    HDIsolution = solve(HDI_prob) # solve for the two roots

    return HDIsolution.u
end