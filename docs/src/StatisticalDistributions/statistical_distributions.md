# Statistical Distributions

The mathematical definitions of several common (in Oceanography) probability density functions and their cumulative distribution functions are defined here for general use.

## General Application
In addition, some generally-useful functions related to statistical distributions are provided.

```@docs; canonical=false
find_pdf_mode
find_HDI
transform_pdf
```

#### Normal Functions
```@docs; canonical=false
gaussian
normal
∫Normaldx
NormalCDF
```

#### Log-Normal Functions
```@docs; canonical=false
lognormal
LogNormalCDF
```

### Skew-Normal Functions
```@docs; canonical=false
skewgaussian
skewnormal
SkewNormalCDF
OwensTFunction
skew_δ
skew_mean
SkewNormalμ2λ
skew_σ
SkewNormalσ2ω
skew_skew
SkewNormalSkewness2Posδ
SkewNormalSkewness2Posξ
skew_mode
```

### Log-Skew-Normal Functions
```@docs; canonical=false
logskewnormal
LogSkewNormalCDF
```