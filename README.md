# Latent binomial correlation

R code for estimating latent covariance and correlation from multivariate
binomial counts with heterogeneous denominators, including simulation and
Australian Electoral Commission (AEC) analyses.

## Requirements

R 4.4 or later.

```r
install.packages("renv")
renv::restore()
```

Core functions are under `R/`; numbered analysis scripts are under `scripts/`
and should be run from the repository root.

The AEC acquisition scripts retrieve public official data directly from the
Commission website.

## Licence

The author-created software is released under the MIT License; see `LICENSE`.
