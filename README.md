# Latent covariance from multivariate binomial counts

Reproducible R code for the manuscript **“Latent covariance and correlation
from multivariate binomial counts with informative denominators”**.

Manuscript prepared for submission to the *Australian & New Zealand Journal
of Statistics*. The manuscript is not accepted or published.

## Statistical problem

An observed proportion combines a latent probability with binomial sampling
variation. Correlating observed proportions can therefore attenuate the
association between the probabilities, especially when denominators are small
or heterogeneous. The implementation estimates the probability-scale
covariance from one aggregate count per unit-coordinate while allowing the
denominators to be mutually dependent and associated with the probability
vector.

For independent units, the working observation model is

```text
X_ij | (P_i, M_i) ~ Binomial(M_ij, P_ij),
```

conditionally independently across coordinates. The primary estimand is
`cov(P_i)` and its ordinary pairwise correlations. No parametric distribution
is imposed on the joint latent law of probabilities and denominators.

## Software

- R 4.4.1, or a compatible R 4.4 installation;
- packages recorded in `renv.lock`;
- internet access for the optional Australian Electoral Commission (AEC)
  acquisition workflow.

Restore the project library from the repository root:

```r
install.packages("renv")
renv::restore(prompt = FALSE)
```

All public code uses project-relative paths. Raw data and large simulation
objects are created locally and are ignored by Git.

## Minimal estimator example

```r
source("R/constants.R")
source("R/utilities.R")
source("R/validate_counts.R")
source("R/psd_projection.R")
source("R/latent_covariance.R")

X <- matrix(c(7, 11, 4, 8, 13, 15, 6, 10), ncol = 2, byrow = TRUE)
M <- matrix(c(20, 25, 12, 18, 30, 32, 16, 22), ncol = 2, byrow = TRUE)

fit <- fit_latent_binomial_cov(X, M)
fit$Sigma_raw       # exactly unbiased latent-covariance estimator
fit$Sigma_psd       # positive-semidefinite projection
fit$Sigma_stable    # n^(-2) eigenvalue floor used for correlations
fit$R_latent        # estimated latent correlation matrix
```

## Quick reproduction

Quick reproduction reads the compact canonical CSV summaries and regenerates
the four main tables, three grayscale vector figures, and supporting tables.
It does not rerun Monte Carlo experiments.

```r
source("reproduce/quick.R")
```

Outputs are written to `reproduce/output/`. On an ordinary laptop this should
take less than one minute.

## Full reproduction

Run the unit tests first:

```r
source("scripts/01_run_unit_tests.R")
```

The individual full workflows are:

```r
# Pairwise simulations and focused interval bootstrap
source("scripts/16_generate_phase3_rng_streams.R")
source("scripts/04_pairwise_full.R")
source("scripts/05_bootstrap_interval_study.R")

# Five-dimensional matrix simulation
source("scripts/18_generate_phase4_rng_streams.R")
source("scripts/07_matrix_full.R")

# AEC acquisition, matching, primary analysis and controlled thinning
source("scripts/08_download_aec.R")
source("scripts/09_prepare_aec.R")
source("scripts/21_adjudicate_aec_outliers.R")
source("scripts/22_freeze_aec_matches.R")
source("scripts/10_analyse_aec.R")
source("scripts/11_controlled_thinning.R")

# Electoral-division enrichment and cluster bootstrap
source("scripts/24_enrich_aec_divisions.R")
source("scripts/25_aec_division_cluster_bootstrap.R")
```

Alternatively, after reading the scripts and acquiring the AEC inputs, set the
explicit opt-in and use the orchestrator:

```r
Sys.setenv(FULL_REPRODUCTION = "YES")
source("reproduce/full.R")
```

Runtimes depend on hardware and download speed. In the validated four-worker
run, the 120,000 pairwise replicates took about six minutes, the 15,000 matrix
replicates about 1.5 minutes, the four focused-bootstrap cells about 17
minutes in total, and the two 4,999-replicate division bootstraps about 1.5
minutes. AEC download and parsing time is network-dependent. Full output can
require several gigabytes; quick reproduction is recommended for reviewers.

## AEC data

The 18 raw AEC CSV files are not redistributed. See
[`data/README.md`](data/README.md) for the official 2019 and 2022 sources,
required files, commands, expected directory structure and provenance.

## Directory structure

```text
R/                  estimators, inference, simulation and AEC functions
scripts/            validated reproduction entry points
config/             locked scenarios, tolerances and master seeds
data/               acquisition documentation; raw files are local only
results/summaries/  compact canonical manuscript summaries
tests/               unit and reproducibility tests
reproduce/           quick and full reader-facing wrappers
```

## Citation

Use the metadata in `CITATION.cff`. The manuscript citation should be updated
with journal publication details if they become available.

## Licence

Author-created software is released under the MIT License; see `LICENSE`.
The licence does not apply to third-party AEC data. Data users must follow the
terms of the official source.
