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
| `"microbiome"` | taxon abundances from shotgun or amplicon sequencing, as raw counts or relative abundance | rare taxa pooled into "Other"; zeros replaced; centred log-ratio (CLR) |
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
- Microbiome: see the next section.
- Repeated measures: samples are assumed independent. Use one sample per
  subject, or regress out subject before the pipeline.

## Microbiome data (`type = "microbiome"`)

Abundances from sequencing are compositional: each sample's values sum to a
fixed total, so correlating log abundances directly creates edges that are not
there. The `"microbiome"` type therefore uses the centred log-ratio (CLR): the
log of each taxon relative to the sample's geometric mean.

```r
res <- run_basil_pipeline(
  data  = list(microbiome = taxa, bile_acids = ba, scfa = scfa),
  types = c(microbiome = "microbiome", bile_acids = "concentration",
            scfa = "concentration"))
```

- Input is a samples x taxa table of raw counts or relative abundances
  (fractions or percent), with 0 for undetected taxa. Both give the same CLR
  values apart from zero replacement: 0.5 for counts, half the sample's
  smallest non-zero value for relative abundance.
- A CLR matrix you computed yourself goes in with `type = "none"`, which uses
  the values as they are (no filtering of rare taxa, no further transform):
  `types = c(microbiome = "none", ...)`. Filter rare taxa before computing the
  CLR; any base of logarithm is fine.
- Use one taxonomic level (species or genus). Mixing levels puts a genus and
  its own species in the same table, and they correlate trivially.
- Taxa with zeros in more than `max_censored` of samples (default 0.5) are not
  analysed, but are pooled into an "Other" bin so the CLR is still taken over
  the whole community. Change it with
  `prep_args = list(microbiome = list(max_censored = 0.2))`.
- Do not use `"counts"` or `"concentration"` for taxa. In a simulation with no
  compositional correction (`"counts"`), about 90% of truly absent taxon-taxon
  edges were selected; with `"microbiome"` about 3%.
- Zeros still cost accuracy. In the same simulation, with about 18% zeros among
  the kept taxa, absent edges stayed at about 3-4% selected, but the 95%
  intervals covered true non-zero edges only about 70% (taxon-taxon) to 78%
  (taxon-metabolite) of the time. Trust which edges are present more than
  their exact size. A stricter filter (`max_censored = 0.2`) helped only
  slightly in that simulation.
- Detection depends on sequencing depth. Drop shallow samples, or add log depth
  to `covariates`.
- If absolute abundances are available (qPCR, spike-in, flow cytometry), scale
  to absolute load and use `"concentration"`; the compositional problem then
  does not arise.
- Taxonomy (family, phylum) and function (SCFA producers, bile-acid
  transformers) make natural feature sets.

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

## Plots

With `out_dir` set, three plots are written:

- `edge_forest.pdf`: the strongest selected edges with their credible intervals.
- `network.pdf`: the 300 strongest selected edges; nodes coloured by assay.
- `edge_heatmap.pdf`: all features against all features, non-selected edges
  blank. The margins carry one track for the assay and one track per feature
  set, marking the set's members. Sets can overlap, so each has its own track.
  By default the 8 largest sets are shown; choose others with
  `heatmap_sets = c("primary_BA", "IL6_signaling")`, or `character(0)` for none.
  Extra tracks of your own are added with `annotation` (next section).

### Your own annotation tracks

`annotation` adds further tracks to the heatmap margins. It is a named list
with one entry per track; each entry is a vector of values named by feature.
Values can be categories or numbers. Features not named are left blank, so a
track may cover one assay only.

```r
ann <- list(
  BA_class    = c(CA = "primary", CDCA = "primary", DCA = "secondary", LCA = "secondary"),
  conjugation = c(GCA = "glycine", TCA = "taurine", CA = "free"),
  phylum      = c(Bacteroides_fragilis = "Bacteroidota", Blautia_obeum = "Bacillota")
)
cols <- list(BA_class = c(primary = "#E69F00", secondary = "#7B3294"))   # optional

res <- run_basil_pipeline(data, types, sets = my_sets, out_dir = "results",
                          annotation = ann, annotation_colors = cols)
```

`annotation_colors` is optional and can cover some tracks only; the rest get
default colours. For a categorical track give one colour per category, named
by category. It also accepts `assay` and feature-set names, to change their
default colours, e.g. `list(assay = c(olink = "steelblue", rna = "tomato"))`
(all assays must then be given a colour).

Annotations describe features, not samples, because both axes of the heatmap
are features. Sample-level variables such as sex or batch cannot be shown
here; they go in `covariates`.

The plots can be redrawn from a result with `plot_edge_forest()`,
`plot_basil_network()` and `plot_edge_heatmap()`, e.g.

```r
h <- plot_edge_heatmap(res$matrices, res$selected, res$prep$features,
                       C = res$C, show_sets = c("primary_BA", "secondary_BA"),
                       annotation = ann, annotation_colors = cols)
grid::grid.newpage(); grid::grid.draw(h$gtable)
```

## Settings guide

Every setting of `run_basil_pipeline()` has a default, so a first run needs
only `data` and `types`. The settings fall into three groups: how the model is
fitted, what an edge is and how it is summarised, and which edges are reported.
Only the first group changes the fit; the other two can be changed afterwards
without refitting, by calling `basil_edges()` and `select_edges()` on the result.

### Settings most likely to matter

| Setting | Default | What it controls | When to change it |
|---|---|---|---|
| `k` | chosen automatically | Number of latent factors, i.e. how many shared patterns of variation the model may use to explain correlations. Too few and real edges are missed or blurred together; too many and noise is modelled as signal. | Leave automatic first. The log prints the chosen value. Then rerun with `k - 1` and `k + 1` and check that the edges you care about survive (see below). |
| `edge_type` | `"correlation"` | Whether an edge is a marginal correlation or a partial correlation (association after accounting for all other features). | Use `"partial"` if you want glasso-like "direct" edges, and treat those intervals as exploratory. |
| `rule` | `"ci"` | How edges are picked for the reported network. `"ci"`: the credible interval excludes zero. `"bfdr"`: the expected share of reported edges that are negligible is at most `alpha`. | Use `"bfdr"` when you report many edges and want a statement about the list as a whole, not edge by edge. |
| `min_abs` | `0.1` | Smallest posterior-mean edge worth reporting. Applied with either rule. | Raise (0.2-0.3) for a sparser, more readable network. Set to 0 to keep every edge that passes the rule. |
| `covariates` | none | Variables removed from every feature before fitting. | Add batch, age, sex or anything else that could create correlations you do not want as edges. |

### Model fitting

| Setting | Default | What it controls | When to change it |
|---|---|---|---|
| `k_max` | components explaining 80% of variance | Largest `k` the automatic search tries. Ignored when `k` is given. | Lower it if the search is slow or picks a value at the top of the range. |
| `n_MC` | `1000` | Number of posterior draws. It does not change the fit, only how precisely interval ends and probabilities are computed. With 1000 draws each end of a 95% interval rests on 25 draws. | 2000-5000 for final results. Run time and memory grow in proportion. |
| `equalize_noise` | `TRUE` | Rescales features once so their unexplained variance is equal, which is what BASIL assumes. See "Calibration and limits". | Set `FALSE` only to reproduce plain BASIL on z-scores. |
| `inverse_normal` | `FALSE` | Replaces each feature's values by normal scores of their ranks, so correlations behave like Spearman correlations. | Turn on if features remain skewed or have outliers after the log transform. |
| `seed` | `1` | Random seed for the posterior draws. | Results are reproducible with the same seed; changing it shows how much Monte Carlo noise there is. |

### Edge summaries

| Setting | Default | What it controls | When to change it |
|---|---|---|---|
| `level` | `0.95` | Width of the credible interval (`lower`, `upper`), and therefore how strict `rule = "ci"` is. | `0.99` for a stricter network, `0.90` for a more liberal one. |
| `delta` | `0.1` | The size below which an edge counts as negligible. Used only for the `prob_gt_delta` and `bfdr` columns, so it matters only with `rule = "bfdr"`. | Set it to the smallest correlation you would consider biologically interesting. |
| `alpha` | `0.05` | Target Bayesian false discovery rate for `rule = "bfdr"`. | `0.01` for a stricter list, `0.10` for a more liberal one. |

### Data preparation

These are per assay and go in `prep_args`, e.g.
`prep_args = list(bile_acids = list(max_censored = 0.3))`.

| Setting | Default | Applies to | What it controls |
|---|---|---|---|
| `max_censored` | `0.5` | all types | Features with more than this fraction below LOQ / below LOD / missing are dropped. Lower it (0.2-0.3) if heavily imputed features worry you: imputed values are identical across samples and distort correlations. |
| `lod` | none | `"concentration"`, `"npx"` | Limit of detection, one number or a named vector by feature. For concentrations, values below it are replaced by `lod/2`. For NPX it is only used to count and filter. |
| `na_is_censored` | `TRUE` | `"concentration"` | `TRUE`: a missing value means "below LOQ". `FALSE`: it means "not measured" and is filled with the feature median. |
| `min_count`, `min_prop` | `10`, `0.2` | `"counts"` | A gene is kept if it has at least `min_count` reads in at least `min_prop` of samples. |

### Fixed by the pipeline

- BASIL's prior settings for the noise variance (`v0 = 1`, `sigma_sq0 = 1`) are
  left at the package defaults. They count as one prior observation against
  n x p data values, so they have no practical effect.
- The shrinkage strengths inside and outside the feature sets (`tau_gamma`,
  `tau_psi`) and the coverage correction (`rho`) are estimated from the data by
  BASIL. The log prints them. A ratio `tau_gamma/tau_psi` above 1 means the
  feature sets explain more of the shared variation than an arbitrary grouping
  would.
- Feature sets with fewer than 2 matched features are dropped, and one
  indicator per assay is added. Call `basil_membership()` directly to change
  either.

### Checking sensitivity to `k`

```r
fits <- lapply(c(4, 5, 6), function(k)
  run_basil_pipeline(data, types, sets = my_sets, k = k, verbose = FALSE))
keys <- lapply(fits, function(r) paste(r$selected$feature1, r$selected$feature2))
length(Reduce(intersect, keys))   # edges selected at every k
sapply(keys, length)              # edges selected at each k
```

Edges that appear at every `k` are the ones to trust.

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
