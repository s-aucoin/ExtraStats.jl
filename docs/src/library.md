# Library

All functions exported by this package are listed below.

## Index
```@index
```


## Functions

### Micellaneous
```@docs
nanmean
nanvar
nanstd
moving_average
bincounts
```

### Least Squares Fitting
```@docs
fittest
OLS
GLS
sqcovar
```

#### Fit Evaluation
```@docs
CalculateRMSE
r²
χ²Test
KuiperTest
```

#### Linear Fit Uncertainty
```@docs
linear_uncertainty
SE_beta
SE_alpha
uncertainty_lines
```


### Statistical Distributions
```@docs
find_pdf_mode
find_HDI
transform_pdf
```

#### Normal Functions
```@docs
gaussian
normal
∫Normaldx
NormalCDF
```

#### Log-Normal Functions
```@docs
lognormal
LogNormalCDF
```

### Skew-Normal Functions
```@docs
skewgaussian
skewnormal
SkewNormalCDF
OwensTFunction
SkewNormalCDF
skew_δ
skew_mean
SkewNormalμ2λ,
skew_σ
SkewNormalσ2ω
skew_skew
SkewNormalSkewness2Posδ
SkewNormalSkewness2Posξ
skew_mode
```

### Log-Skew-Normal Functions
```@docs
logskewnormal
LogSkewNormalCDF
```