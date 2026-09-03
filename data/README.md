# Australian Electoral Commission data

The application uses public official House of Representatives data from the
Australian Electoral Commission (AEC).

Official sources:

- 2019: <https://results.aec.gov.au/24310/Website/HouseDownloadsMenu-24310-Csv.htm>
- 2022: <https://results.aec.gov.au/27966/website/HouseDownloadsMenu-27966-Csv.htm>

## Required files

For each election, the workflow obtains eight state/territory
first-preference-by-polling-place CSV files and one national polling-place
metadata CSV file: 18 files in total. Exact file roles and URLs are listed in
`config/aec_manifest.csv`.

## Commands

From the repository root, after restoring `renv.lock`, run:

```r
source("scripts/08_download_aec.R")
source("scripts/09_prepare_aec.R")
source("scripts/21_adjudicate_aec_outliers.R")
source("scripts/22_freeze_aec_matches.R")
source("scripts/10_analyse_aec.R")
source("scripts/11_controlled_thinning.R")
source("scripts/24_enrich_aec_divisions.R")
source("scripts/25_aec_division_cluster_bootstrap.R")
```

The acquisition script downloads into:

```text
data/raw/aec/2019/
data/raw/aec/2022/
```

Derived local files are written under `results/derived/`, and scientific
summaries under `results/summaries/`. Raw and large derived files are ignored
by Git.

## Provenance and construction

The scripts record the requested URL, official basename, retrieval metadata
and cryptographic hash. Preparation detects and validates source schemas,
retains the official ordinary vote type, sums formal ordinary first-preference
votes over all candidates for the denominator, and sums Australian Labor Party
rows for the numerator. Informal votes are excluded. Matching is exact by state
and polling-place identifier where possible; the restricted secondary route
requires a unique normalised name within state and coordinate support within
two kilometres. No unconstrained fuzzy matching is used.

This repository does not redistribute the 18 source CSV files. Users should
review the AEC source information and applicable terms before downloading.
