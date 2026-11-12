# BASIL network pipeline (metabolomics / Olink / targeted RNA)

One call from raw assay tables to a network whose edges carry posterior
credible intervals.

```r
# once: install BASIL from this checkout (includes a fix for single-set C)
devtools::install(".")          # also needs: matrixStats, ggraph, pheatmap

source("pipeline/basil_pipeline.R")

res <- run_basil_pipeline(
  data  = list(bile_acids = ba, scfa = scfa, olink = npx, rna = counts),
  types = c(bile_acids = "concentration", scfa = "concentration",
            olink = "npx", rna = "normalized"),   # FOLI-seq: "none" if already log-scaled
  sets       = my_sets,      # optional named list of feature sets
  covariates = meta,         # optional data.frame: batch, age, sex, ...
  out_dir    = "results"
)
res$selected                 # edges with mean, lower, upper, probabilities
```

Every assay is a **samples x features** numeric matrix with sample IDs as
rownames. Assays are matched on sample ID; only samples present in all of them
are used. A single assay works too. `Rscript pipeline/example_simulated.R` runs
a complete example.

## What goes in and what comes out

BASIL does **not** take a p x p matrix. It takes

- `Y`: the n x p data matrix (the pipeline builds it), and
- `C`: a p x q 0/1 matrix saying which feature belongs to which set (pathway,
  chemical class, ...). The pipeline builds it from `sets`.

It fits a factor model, `Cov(Y) = Lambda Lambda' + sigma^2 I`, and returns
posterior draws of `Lambda` and `sigma^2`. The p x p matrix is the *output*:
each draw gives a covariance matrix, hence a posterior distribution for every
pairwise edge.

| `edge_type` | Edge means | Status |
|---|---|---|
| `"correlation"` (default) | marginal correlation, includes indirect association | what the BASIL paper and its coverage correction address |
| `"partial"` | partial correlation from the implied precision matrix (comparable to glasso edges) | exploratory: derived per draw by the pipeline; BASIL's coverage correction was not designed for it, and the precision matrix of a factor model is dense, not sparse |

## Data preparation (done by `types`)

| type | Use for | What happens |
|---|---|---|
| `"concentration"` | bile acids, SCFAs (absolute or relative concentrations) | zeros / below-`lod` / NA treated as below LOQ and set to half the feature's smallest observed value (or `lod/2`); log2 |
| `"npx"` | Olink NPX | no transform (NPX is log2); values below LOD are kept; sporadic NAs median-imputed |
| `"counts"` | raw read counts | low-count filter (>= 10 reads in >= 20% of samples), log2-CPM |
| `"normalized"` | expression already normalised but on the linear scale (normalised counts, CPM/TPM) | log2(x + 1) |
| `"none"` | anything already normalised and log-scaled | missing-value handling only |

All types drop features with > 50% censored/missing (`max_censored`) or no
variance. Per-assay options go in `prep_args`, e.g.
`prep_args = list(bile_acids = list(lod = lod_vector, na_is_censored = FALSE))`.
Then, for all features together: optional covariate regression, optional
rank-based inverse normal transform (`inverse_normal = TRUE`, for skewed or
outlier-prone data), z-scoring.

Do upstream, before the pipeline:

- Olink: plate/batch bridging and sample QC (e.g. `OlinkAnalyze`).
  `olink_long_to_wide()` converts the long NPX export to a matrix.
- Metabolomics: normalise stool concentrations to sample weight if not done.
- Expression that is already normalised (e.g. FOLI-seq): use `"normalized"` if
  values are on the linear scale (non-negative, right-skewed, maximum in the
  thousands) and `"none"` if they are already log-transformed. `"counts"` is
  only for raw reads; on a small targeted panel its CPM step uses the panel
  total as library size, which is poor when a few genes dominate.
- Repeated measures: samples are assumed independent. Use one sample per
  subject, or regress out subject before the pipeline.

## Feature sets (`sets`)

A named list of character vectors; a set may span assays.

```r
my_sets <- list(
  primary_BA    = c("CA", "CDCA", "GCA", "TCA", ...),
  secondary_BA  = c("DCA", "LCA", "UDCA", ...),
  conjugated_BA = c("GCA", "TCA", "GDCA", ...),
  IL6_signaling = c("IL6", "IL6R", "SOCS3", ...)   # proteins and genes together
)
```

With several assays, an indicator per assay is added automatically. Without
`sets` the fit is a plain Bayesian factor model: still valid, but BASIL's
pathway-informed shrinkage is not used.

## Reading the edge table

| column | meaning |
|---|---|
| `mean`, `lower`, `upper` | posterior mean and 95% credible interval of the edge |
| `prob_direction` | posterior probability that the edge has the sign of its mean |
| `prob_gt_delta` | posterior probability that \|edge\| > `delta` (default 0.1) |
| `bfdr` | Bayesian FDR if all edges at least this certain are reported |
| `ci_excludes_zero` | interval does not contain 0 |

Selection (`rule`): `"ci"` keeps edges whose interval excludes zero and
\|mean\| >= `min_abs` (default 0.1); `"bfdr"` keeps edges with `bfdr <= alpha`.
`edges_all.csv` has every pair, so other rules can be applied later with
`select_edges()`.

## Calibration and limits

Checked on simulated factor-model data with unequal noise across features and
30% pure-noise features (n = 60-300, p = 150-300):

- 95% intervals cover true non-zero correlations about 91-92% of the time and
  true zero correlations about 97%.
- About 1-3% of truly absent edges are selected with the default rule.
- Without `equalize_noise` (plain z-scores) coverage of non-zero edges drops
  to about 80% and strong edges are biased toward zero. `equalize_noise = TRUE`
  is a pipeline addition: one rescaling of features to equal estimated residual
  variance, because BASIL assumes a single `sigma^2`.

Limits:

- Intervals are per edge; nothing adjusts for the number of edges except the
  `bfdr` rule.
- The model needs real low-rank structure. If JIC picks k = 1 the pipeline
  fits k = 2 (BASIL cannot fit k = 1) and reports it; on pure-noise data that
  produced spurious weak edges (about 7% of pairs, all near \|r\| = 0.1-0.2).
- Calibration was checked on Gaussian simulations, not on real assay data.
- A large assay (e.g. Olink Explore, ~3000 proteins) dominates the factors when
  combined with a 20-feature panel; edges within the small panel are then
  estimated mainly through what it shares with the large one.
- Covariate regression does not reduce the sample size used by the model, so
  intervals are slightly narrow when many covariates are removed.
