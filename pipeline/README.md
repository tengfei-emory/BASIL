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

## Method

This section describes what the code does, step by step. Steps 2-4 are the
BASIL package; steps 1, 5 and 6 are added by the pipeline.

### 1. Preparing the data

Each assay is transformed to a roughly Gaussian scale (see the next section),
optionally residualised on covariates by least squares, and every feature is
centred and scaled to unit variance. The assays are joined into one matrix `Y`
with n samples (those present in all assays) and p features.

### 2. The model

BASIL is a Bayesian latent factor model. Each sample is explained by `k`
unobserved factors plus independent noise:

```
y_i = Lambda eta_i + e_i,   eta_i ~ N(0, I_k),   e_i ~ N(0, sigma^2 I_p)
Cov(y_i) = Lambda Lambda' + sigma^2 I_p
```

`Lambda` (p x k) holds the loadings: how strongly each feature follows each
factor. Two features are correlated when they load on the same factors. With
k much smaller than p, the whole p x p covariance matrix is described by
p x k numbers, which is what makes it estimable when p is larger than n.

The feature sets enter through the prior on `Lambda`. The loadings are split
into a part that is constant within sets (the part lying in the column space
of the membership matrix `C`) and a remainder. Each part has its own Gaussian
prior, with variance `tau_gamma sigma^2` and `tau_psi sigma^2`. Both variances
are estimated from the data (empirical Bayes). If the sets explain much of the
shared variation, `tau_gamma` comes out larger than `tau_psi`, and the
set-driven part is shrunk less than the rest. Uninformative sets get a small
`tau_gamma`, so their part is shrunk strongly.

### 3. Fitting: no MCMC

BASIL avoids Markov chain Monte Carlo. It works in two moves:

1. **Estimate the factors.** The factor scores are taken from the leading `k`
   singular vectors of `Y`. With many features this estimate is accurate,
   because every feature carries information about the same few factors.
   `k` is chosen by an information criterion (JIC) unless you fix it.
2. **Posterior of the loadings given the factors.** With the factors treated
   as known, the model is a set of linear regressions of each feature on the
   factors, and the posterior is available in closed form:
   - `sigma^2` has an inverse-gamma posterior;
   - given `sigma^2`, `Lambda` is Gaussian, centred at a shrunken regression
     estimate (shrinkage factor `n / (n + 1/tau)`, separately for the two
     parts), with standard deviation `sigma / sqrt(n + 1/tau)`.

Posterior draws are therefore independent draws from known distributions.
There is no burn-in and no convergence to check, and a fit takes seconds.

### 4. Coverage correction

Treating the estimated factors as known ignores one source of uncertainty, so
the closed-form posterior is too narrow. BASIL corrects this by multiplying
the posterior standard deviation of the loadings by a factor `rho`, which is
at least 1 and is computed from the fitted loadings and noise variance. The
purpose is that credible intervals for the entries of the covariance matrix
cover the true value at the stated rate. The log prints `rho`; it was about
1.2-1.3 in our simulated examples.

### 5. From loadings to edges with uncertainty

The pipeline takes `n_MC` posterior draws (default 1000) of `(Lambda,
sigma^2)`. For each draw `t` it forms the implied matrix

```
Sigma_t = Lambda_t Lambda_t' + sigma_t^2 I
```

and converts it to correlations (or, with `edge_type = "partial"`, inverts it
and converts to partial correlations). Every feature pair `(i, j)` thus gets
`n_MC` values of its edge: a sample from the posterior distribution of that
edge. From this sample the pipeline reports

| quantity | how it is computed |
|---|---|
| `mean` | average over draws |
| `lower`, `upper` | the 2.5% and 97.5% quantiles of the draws (for `level = 0.95`) |
| `prob_direction` | fraction of draws with the same sign as the mean |
| `prob_gt_delta` | fraction of draws with \|edge\| > `delta` |

This is the difference from glasso. Glasso returns one sparse estimate, and
an edge is either in or out with no statement of confidence. Here every edge
has a full posterior distribution, so a weak but precisely estimated edge can
be told apart from a strong but uncertain one.

### 6. Selecting edges

- `rule = "ci"`: report an edge if its credible interval excludes zero and
  \|mean\| >= `min_abs`. Each interval is a per-edge statement; nothing
  accounts for the number of edges examined.
- `rule = "bfdr"`: Bayesian false discovery rate. Edges are ranked by
  `prob_gt_delta`. For an edge, `1 - prob_gt_delta` is the posterior
  probability that it is *not* larger than `delta`. The `bfdr` of an edge is
  the average of this probability over the edge and all edges ranked above
  it, which is the expected fraction of false reports if the list is cut
  there. Edges with `bfdr <= alpha` are reported. This rule does account for
  the number of edges.

### What the uncertainty does and does not include

Included: sampling variability of the loadings and of `sigma^2`, with the
correction for estimated factors.

Not included: the choice of `k`, the estimated prior variances, the data
preparation (imputation, covariate adjustment, CLR), and the `equalize_noise`
rescaling. All of these are treated as fixed. The intervals also rest on the
model: Gaussian data, independent samples, a common noise variance, and
low-rank structure. "Calibration and limits" below reports how the intervals
behave in simulations.

The model, the fitting procedure and the coverage correction are those of
BASIL. For the theory and the simulation studies, see the paper:

> Mauri L, Stolf F, Herring AH, Miller C, Dunson DB (2026). Pathway-based
> Bayesian factor models for 'omics data. arXiv:2601.13419.
> <https://arxiv.org/abs/2601.13419>

Please cite it when you report results from this pipeline. The package is at
<https://github.com/federicastolf/BASIL>.

## Data preparation (done by `types`)

| type | Use for | What happens |
|---|---|---|
| `"concentration"` | bile acids, SCFAs (absolute or relative concentrations) | zeros / below-`lod` / NA treated as below LOQ and set to half the feature's smallest observed value (or `lod/2`); log2 |
| `"npx"` | Olink NPX | no transform (NPX is log2); values below LOD are kept; sporadic NAs median-imputed |
| `"counts"` | raw read counts | low-count filter (>= 10 reads in >= 20% of samples), log2-CPM; then the 2000 most variable genes |
| `"normalized"` | expression already normalised but on the linear scale (normalised counts, CPM/TPM) | log2(x + 1); then the 2000 most variable genes |
| `"microbiome"` | taxon abundances from shotgun or amplicon sequencing, as raw counts or relative abundance | rare taxa pooled into "Other"; zeros replaced; centred log-ratio (CLR) |
| `"none"` | anything already normalised and log-scaled | missing-value handling only |

All types drop features with > 50% censored/missing (`max_censored`) or no
variance. Per-assay options go in `prep_args`, e.g.
`prep_args = list(bile_acids = list(lod = lod_vector, na_is_censored = FALSE))`.

**Highly variable genes.** A whole transcriptome should not go into the
network: most genes barely vary between samples, and they add noise and
millions of uninformative pairs. For `"counts"` and `"normalized"` the
pipeline therefore keeps the `n_hvg` most variable genes (default 2000, a
common choice). Genes are ranked on the log scale by their variance relative
to the mean-variance trend (a loess fit of log variance on mean expression),
the same idea as the variable-gene selection in Seurat and scran: a gene
counts as variable if it varies more than other genes of similar expression
level, so genes are not picked just because low expression makes them noisy.
A panel with fewer genes than `n_hvg` is left as it is. Change the number
with `prep_args = list(rna = list(n_hvg = 1000))`; `n_hvg = Inf` keeps all
genes. The other types have no limit by default, but accept `n_hvg` too
(useful for log-scale expression entered as `"none"`).
`feature_prep_log.csv` records which genes were dropped and why.

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

A set that contains exactly the features of one assay (or exactly the members
of another set) carries no extra information for the model, so it is left out
of the fit. The log says so, and the set still gets its track in the heatmap.

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
  With several assays, features are grouped by assay with a gap between
  assays, so each off-diagonal block is one pair of assays; clustering is
  done within each assay. `heatmap_split = FALSE` clusters all features
  together instead.
  Feature names are printed for up to 200 features; above 60 the font
  shrinks and the page grows, so zoom in. Force them on or off with
  `heatmap_names = TRUE` / `FALSE`.

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
| `n_hvg` | `2000` for `"counts"` and `"normalized"`, no limit otherwise | all types | Keep only this many of the most variable features (highly variable genes). Fewer genes give a faster fit and a more readable network; `Inf` keeps all. Genes named in `sets` but filtered out are simply absent from the analysis. |

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
