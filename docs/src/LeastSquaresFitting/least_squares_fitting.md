# Least Squares Fitting

These functions are specific to least squares fitting applications.

## General Application
Least squares fitting using ordinary or generalized least squares is supported.

```@docs; canonical=false
fittest
OLS
GLS
sqcovar
```


## Fit Evaluation
Several functions for evaluating the goodness-of-fit of a given fit are provided.

```@docs; canonical=false
CalculateRMSE
r²
χ²Test
KuiperTest
```


## Linear Fit Uncertainty
Specifically for linear fits, uncertainty on the slope and intercept can easily be estimated.

```@docs; canonical=false
linear_uncertainty
SE_beta
SE_alpha
uncertainty_lines
```