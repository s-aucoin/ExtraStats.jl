module LinearFitUncertainty

using Statistics
using Distributions

export linear_uncertainty, SE_beta, SE_alpha, uncertainty_lines


"""
    linear_uncertainty(xi, yi, yhat, γ)

Calculate the confidence intervals at the `γ` level on the fit parameters of a linear regression on indpendent observations `xi`, dependent observations `yi`, and predicted dependent values `yhat`.
"""
function linear_uncertainty(xi, yi, yhat, γ)

    n = size(xi, 1)
    dof = n - 2

    t_val = quantile.(TDist(dof), 1 - γ/2)

        beta_e = SE_beta(xi, yi, yhat) .* t_val
    alpha_e = SE_alpha(xi, yi, yhat) .* t_val

    return (; beta_e, alpha_e)
end



"""
    SE_beta(xi, yi, yhat)

Calculate the standard error of the slope of a linear regression with independent observations `xi`, dependent observations `yi`, and predicted dependent values `yhat`.
"""
function SE_beta(xi, yi, yhat)
    n = size(xi)[1]
    xbar = mean(xi)
    return sqrt(1/(n-2) .* sum((yi .- yhat).^2) ./ (sum((xi .- xbar).^2)))
end


"""
    SE_alpha(xi, yi, yhat)

Calculate the standard error of the intercept of a linear regression with independent observations `xi`, dependent observations `yi`, and predicted dependent values `yhat`.
"""
function SE_alpha(xi, yi, yhat)
    n = size(xi)[1]
    s_beta = SE_beta(xi, yi, yhat)
    return s_beta .* sqrt(1/n .* sum(xi.^2))
end



"""
    uncertainty_lines(β_est, α_est, x, t_val, xi, yi, yhat)

Create the upper and lower confidence bound lines of a linear regression.

...
# Arguments
- `β_est`: best fit estimate of the slope.
- `α_est`: best fit estimate of the intercept.
- `x`: independent variable values to calculate bounds at.
- `xi`: independent variable observations.
- `yi`: dependent variable observations.
- `yhat`: predicted values of the dependent variable.
- `γ`: confidence level (i.e. γ-th percentile).
...
"""
function uncertainty_lines(β_est, α_est, x, xi, yi, yhat, γ)
    n = size(xi)[1]
    xbar = mean(xi)

    dof = n - 2
    t_val = quantile.(TDist(dof), 1 - γ/2)

    ξ = t_val .* sqrt.(1/(n-2)*sum((yi .- yhat).^2)*(1/n .+ (x .- xbar).^2 ./ (sum((xi .- xbar).^2))))

    upper_ln = α_est .+ β_est .* x .+ ξ
    lower_ln = α_est .+ β_est .* x .- ξ

    return (; upper_ln, lower_ln)
end


end