# =============================================================================
# BASIL network pipeline for (multi-)omics data
#
#   source("pipeline/basil_pipeline.R")
#   res <- run_basil_pipeline(data = list(...), types = c(...), out_dir = "results")
#
# Steps (each is also a stand-alone function):
#   1. prep_modality()     per-assay cleaning / transformation
#   2. basil_prepare()     match samples, (optionally) adjust covariates, z-score
#   3. basil_membership()  build the feature-set matrix C that BASIL needs
#   4. basil_fit()         choose k, fit BASIL, draw posterior samples
#   5. basil_edges()       per-edge posterior mean, credible interval, probabilities
#   6. plot_*()            forest plot of edges, network, heatmap
#
# Requires: BASIL, matrixStats, ggplot2, igraph, ggraph (the last three are
# BASIL dependencies already).
# =============================================================================

.need <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required. Install it first.", call. = FALSE)
}
.msg <- function(verbose, ...) if (verbose) message(...)

# -----------------------------------------------------------------------------
# 1. Per-assay preparation
# -----------------------------------------------------------------------------

#' Clean and transform one assay
#'
#' @param x samples x features numeric matrix / data.frame, rownames = sample IDs.
#' @param type
#'   "concentration" targeted metabolomics (bile acids, SCFAs). Zeros, values
#'                   below `lod` and (if `na_is_censored`) NAs are treated as
#'                   left-censored and replaced by half the smallest observed
#'                   value of that feature (or lod/2), then log2.
#'   "npx"           Olink NPX. Already log2; no transformation. Sporadic NAs
#'                   are median-imputed.
#'   "counts"        raw sequencing counts. Low-count filter, then log2-CPM.
#'   "normalized"    already normalised but still on the linear scale
#'                   (normalised counts, CPM/TPM, ...): log2(x + 1) only.
#'   "none"          already normalised and on a log-like scale.
#' @param lod optional limit of detection: one number or a named vector by feature.
#' @param max_censored drop features with more than this fraction censored /
#'   below LOD / missing (default 0.5).
#' @param na_is_censored for "concentration": are NAs "below LOQ" (TRUE) or
#'   "not measured" (FALSE -> median-imputed on the log scale)?
#' @param min_count,min_prop for "counts": keep features with at least
#'   `min_count` reads in at least `min_prop` of samples.
#' @return list(x = cleaned matrix, log = per-feature record of what was done)
prep_modality <- function(x, type = c("concentration", "npx", "counts", "normalized", "none"),
                          lod = NULL, max_censored = 0.5, na_is_censored = TRUE,
                          min_count = 10, min_prop = 0.2) {
  type <- match.arg(type)
  x <- as.matrix(x)
  if (!is.numeric(x)) stop("Assay matrix must be numeric (samples x features).")
  if (is.null(rownames(x)) || is.null(colnames(x)))
    stop("Assay matrix needs rownames (sample IDs) and colnames (feature names).")
  if (anyDuplicated(rownames(x))) stop("Duplicated sample IDs.")
  if (anyDuplicated(colnames(x)))
    stop("Duplicated feature names: ",
         paste(head(unique(colnames(x)[duplicated(colnames(x))])), collapse = ", "))
  p0 <- ncol(x)
  lod_vec <- rep(NA_real_, p0); names(lod_vec) <- colnames(x)
  if (!is.null(lod)) {
    if (length(lod) == 1 && is.null(names(lod))) lod_vec[] <- lod
    else lod_vec[intersect(names(lod), colnames(x))] <-
        lod[intersect(names(lod), colnames(x))]
  }
  below <- sweep(x, 2, lod_vec, "<"); below[is.na(below)] <- FALSE
  log <- data.frame(feature = colnames(x), frac_missing = colMeans(is.na(x)),
                    frac_censored = 0, kept = TRUE, reason = "",
                    stringsAsFactors = FALSE)

  if (type == "concentration") {
    if (any(x < 0, na.rm = TRUE))
      stop("Negative values in a 'concentration' assay. Is it already log-transformed? Use type = 'none'.")
    cens <- below | (!is.na(x) & x <= 0)
    if (na_is_censored) cens <- cens | is.na(x)
    log$frac_censored <- colMeans(cens)
    for (j in seq_len(p0)) {
      obs <- x[!cens[, j] & !is.na(x[, j]), j]
      fill <- if (!is.na(lod_vec[j])) lod_vec[j] / 2
              else if (length(obs)) min(obs) / 2 else NA_real_
      x[cens[, j], j] <- fill
    }
    x <- log2(x)
  } else if (type == "npx") {
    log$frac_censored <- colMeans(below)   # below LOD: counted, values kept
  } else if (type == "counts") {
    if (anyNA(x)) stop("Count matrix contains NA.")
    if (any(x < 0)) stop("Negative values in a count matrix.")
    lib <- rowSums(x)
    if (any(lib == 0)) stop("Samples with zero total counts: ",
                            paste(rownames(x)[lib == 0], collapse = ", "))
    low <- colMeans(x >= min_count) < min_prop
    log$kept[low] <- FALSE; log$reason[low] <- "low counts"
    x <- log2((x + 0.5) / (lib + 1) * 1e6)
  } else if (type == "normalized") {
    if (any(x < 0, na.rm = TRUE))
      stop("Negative values in a 'normalized' assay: it looks log-transformed already. Use type = 'none'.")
    x <- log2(x + 1)
  }

  drop <- log$kept & pmax(log$frac_censored, log$frac_missing *
                            !(type == "concentration" && na_is_censored)) > max_censored
  log$kept[drop] <- FALSE; log$reason[drop] <- "too many censored/missing"

  # remaining NAs: median on the (log) analysis scale
  for (j in which(log$kept & colSums(is.na(x)) > 0))
    x[is.na(x[, j]), j] <- stats::median(x[, j], na.rm = TRUE)

  v <- apply(x, 2, stats::var)
  flat <- log$kept & (is.na(v) | v < 1e-12)
  log$kept[flat] <- FALSE; log$reason[flat] <- "no variance"

  list(x = x[, log$kept, drop = FALSE], log = log)
}

#' Reshape an Olink long-format NPX table to samples x assays
#'
#' @param npx data.frame with columns SampleID, NPX and an assay identifier.
#' @param id column used as feature name. "Assay" is readable, but the same
#'   protein can sit on several panels (e.g. IL6, TNF, CXCL8 in Explore); those
#'   become "Assay|Panel". Use id = "OlinkID" to keep everything distinct.
#' @param drop_qc_warning drop rows whose QC_Warning is not "PASS".
olink_long_to_wide <- function(npx, id = "Assay", drop_qc_warning = TRUE) {
  npx <- as.data.frame(npx)
  stopifnot(all(c("SampleID", "NPX", id) %in% names(npx)))
  if (drop_qc_warning && "QC_Warning" %in% names(npx))
    npx <- npx[toupper(npx$QC_Warning) == "PASS", , drop = FALSE]
  feat <- as.character(npx[[id]])
  if ("Panel" %in% names(npx)) {
    multi <- tapply(as.character(npx$Panel), feat, function(z) length(unique(z)) > 1)
    is_multi <- feat %in% names(multi)[multi]
    feat[is_multi] <- paste(feat[is_multi], npx$Panel[is_multi], sep = "|")
  }
  w <- tapply(npx$NPX, list(as.character(npx$SampleID), feat), mean, na.rm = TRUE)
  w[is.nan(w)] <- NA
  as.matrix(w)
}

# -----------------------------------------------------------------------------
# 2. Combine assays into the BASIL input matrix
# -----------------------------------------------------------------------------

.inverse_normal <- function(v) stats::qnorm((rank(v) - 0.5) / length(v))

#' Prepare the samples x features matrix Y
#'
#' @param data named list of assays (samples x features, rownames = sample IDs).
#'   One assay is fine. Only samples present in every assay are used.
#' @param types named or positional character vector, one of prep_modality()'s
#'   types per assay.
#' @param prep_args optional named list of extra prep_modality() arguments per
#'   assay, e.g. list(bile_acids = list(lod = 0.01)).
#' @param covariates optional data.frame (rownames = sample IDs) of nuisance
#'   variables (batch, age, sex, ...) regressed out of every feature.
#' @param inverse_normal rank-based inverse normal transform per feature
#'   (robust to outliers / skew; correlations become Spearman-like).
#' @return list(Y, features, prep_log, samples)
basil_prepare <- function(data, types, prep_args = list(), covariates = NULL,
                          inverse_normal = FALSE, verbose = TRUE) {
  if (!is.list(data) || is.data.frame(data) || is.null(names(data)) || any(names(data) == ""))
    stop("`data` must be a *named* list of assays, e.g. list(olink = npx_matrix).")
  if (is.null(names(types))) {
    if (length(types) != length(data)) stop("`types` needs one entry per assay.")
    names(types) <- names(data)
  }
  if (!all(names(data) %in% names(types))) stop("`types` is missing some assays.")

  ids <- Reduce(intersect, lapply(data, rownames))
  if (!is.null(covariates)) {
    covariates <- as.data.frame(covariates)
    ok <- rownames(covariates)[stats::complete.cases(covariates)]
    ids <- intersect(ids, ok)
  }
  if (length(ids) < 10)
    stop("Only ", length(ids), " samples shared across assays",
         if (!is.null(covariates)) " and covariates", ". Check that rownames are sample IDs.")
  .msg(verbose, "Samples shared by all assays: ", length(ids))

  mats <- list(); logs <- list()
  for (m in names(data)) {
    args <- c(list(x = as.matrix(data[[m]])[ids, , drop = FALSE], type = types[[m]]),
              prep_args[[m]])
    pm <- do.call(prep_modality, args)
    pm$log$modality <- m
    logs[[m]] <- pm$log
    mats[[m]] <- pm$x
    .msg(verbose, sprintf("  %-12s (%s): kept %d of %d features", m, types[[m]],
                          ncol(pm$x), nrow(pm$log)))
  }
  features <- do.call(rbind, lapply(names(mats), function(m)
    data.frame(feature = colnames(mats[[m]]), name = colnames(mats[[m]]),
               modality = m, stringsAsFactors = FALSE)))
  dup <- features$feature %in% features$feature[duplicated(features$feature)]
  features$feature[dup] <- paste(features$modality[dup], features$name[dup], sep = ":")
  Y <- do.call(cbind, mats); colnames(Y) <- features$feature

  if (!is.null(covariates)) {
    X <- stats::model.matrix(~ ., data = covariates[ids, , drop = FALSE])
    Y <- stats::lm.fit(X, Y)$residuals
    dimnames(Y) <- list(ids, features$feature)
    .msg(verbose, "Regressed out: ", paste(names(covariates), collapse = ", "))
  }
  if (inverse_normal) Y <- apply(Y, 2, .inverse_normal)
  Y <- scale(Y)
  attr(Y, "scaled:center") <- attr(Y, "scaled:scale") <- NULL
  rownames(Y) <- ids
  if (anyNA(Y)) stop("NA after scaling (constant feature after covariate adjustment?).")
  rownames(features) <- NULL
  list(Y = Y, features = features, prep_log = do.call(rbind, logs), samples = ids)
}

# -----------------------------------------------------------------------------
# 3. Feature-set matrix C
# -----------------------------------------------------------------------------

#' Build BASIL's feature x set membership matrix
#'
#' BASIL shrinks loadings differently inside / outside the column space of C,
#' so C is where prior biology goes (bile-acid classes, protein pathways, gene
#' sets, ...). With `sets = NULL` only assay indicators (or a single all-ones
#' column) are used and BASIL behaves like a plain Bayesian factor model.
#'
#' @param features `features` data.frame from basil_prepare().
#' @param sets named list of character vectors. Members are matched against
#'   `features$feature` and `features$name`. A set may span assays.
#' @param add_modality add one indicator column per assay.
#' @param min_size drop sets with fewer matched members.
basil_membership <- function(features, sets = NULL, add_modality = TRUE,
                             min_size = 2, verbose = TRUE) {
  p <- nrow(features)
  cols <- list()
  if (add_modality || is.null(sets))
    for (m in unique(features$modality))
      cols[[paste0("assay:", m)]] <- as.numeric(features$modality == m)
  if (!is.null(sets)) {
    if (!is.list(sets) || is.null(names(sets))) stop("`sets` must be a named list.")
    n_drop <- 0
    for (s in names(sets)) {
      v <- as.numeric(features$feature %in% sets[[s]] | features$name %in% sets[[s]])
      if (sum(v) >= min_size && sum(v) < p) cols[[s]] <- v else n_drop <- n_drop + 1
    }
    .msg(verbose, "Feature sets used: ", length(sets) - n_drop, " of ", length(sets),
         " (", n_drop, " dropped: fewer than ", min_size, " matched members)")
  }
  C <- do.call(cbind, cols)
  rownames(C) <- features$feature
  C <- C[, !duplicated(t(C)), drop = FALSE]
  q <- qr(C)$rank
  if (q >= p)
    stop("C has rank ", q, " >= number of features (", p, "). Use fewer / broader sets.")
  .msg(verbose, "Membership matrix C: ", p, " features x ", ncol(C), " sets (rank ", q, ")")
  C
}

# -----------------------------------------------------------------------------
# 4. Fit
# -----------------------------------------------------------------------------

#' Fit BASIL and draw posterior samples of the loadings
#'
#' @param Y scaled samples x features matrix (basil_prepare()$Y).
#' @param C membership matrix (basil_membership()).
#' @param k number of latent factors; NULL = choose by JIC.
#' @param k_max largest k tried by JIC; NULL = number of components explaining
#'   80% of variance (as in the BASIL vignette).
#' @param n_MC number of posterior draws.
#' @param equalize_noise BASIL assumes one residual variance for all features.
#'   After z-scoring that only holds if every feature is equally "explained" by
#'   the factors. If TRUE (default), features are rescaled once so that their
#'   estimated residual variance is equal before the final fit. Correlations
#'   are unaffected by the rescaling itself; in simulations it removes the
#'   attenuation of strong edges and brings interval coverage close to nominal.
#'   (Pipeline addition, not part of the BASIL package. FALSE = plain z-scores.)
basil_fit <- function(Y, C, k = NULL, k_max = NULL, n_MC = 1000, seed = 1,
                      equalize_noise = TRUE, verbose = TRUE) {
  .need("BASIL")
  Y <- as.matrix(Y); C <- as.matrix(C); storage.mode(C) <- "double"
  n <- nrow(Y); p <- ncol(Y)
  if (nrow(C) != p) stop("nrow(C) must equal ncol(Y).")
  if (p < 4 || n < 5) stop("Too few features or samples.")
  k_cap <- min(n, p) - 2          # BASIL needs at least 2 residual directions
  jic <- NULL
  if (is.null(k)) {
    if (is.null(k_max)) {
      d2 <- svd(Y, nu = 0, nv = 0)$d^2
      k_max <- which(cumsum(d2) / sum(d2) > 0.8)[1]
    }
    k_max <- max(2, min(k_max, k_cap))
    grDevices::pdf(NULL)           # estimate_latent_dimension() draws a plot
    est <- tryCatch({
      utils::capture.output(e <- BASIL::estimate_latent_dimension(Y, k_max = k_max))
      e
    }, finally = grDevices::dev.off())
    jic <- est$jics; k <- est$k_hat
    .msg(verbose, "JIC chose k = ", k, " (searched 1..", k_max, ")")
  }
  if (k < 2) {                     # k = 1 is not supported by BASIL's matrix code
    .msg(verbose, "Using k = 2 (BASIL does not support k = 1).")
    k <- 2
  }
  if (k > k_cap) stop("k = ", k, " is too large for n = ", n, ", p = ", p, ".")
  if (equalize_noise) {
    s <- svd(Y, nu = k, nv = k)
    res <- Y - s$u %*% (s$d[1:k] * t(s$v))
    Y <- sweep(Y, 2, sqrt(colSums(res^2) / (n - k)), "/")
  }

  set.seed(seed)
  utils::capture.output(
    fit <- BASIL::BASIL_point_estimates(Y = Y, C = C, k = k))
  utils::capture.output(
    ps <- BASIL::compute_posterior_samples_cc(
      Y = Y, Lambda_C = fit$Lambda_C, Lambda_N = fit$Lambda_N,
      tau_gamma = fit$tau_gamma, tau_psi = fit$tau_psi, sigma_sq = fit$sigma_sq,
      P_C = fit$P_C, n_MC = n_MC))
  .msg(verbose, sprintf("BASIL fit: n = %d, p = %d, k = %d, sigma^2 = %.3f, tau_gamma/tau_psi = %.2f, rho = %.2f",
                        n, p, k, fit$sigma_sq, fit$tau_gamma / fit$tau_psi, ps$rho))
  list(fit = fit, Lambda_samples = ps$Lambda_samples,
       sigma_sq_samples = as.numeric(ps$sigma_sq_samples), rho = ps$rho,
       k = k, jic = jic, features = colnames(Y), n = n)
}

# -----------------------------------------------------------------------------
# 5. Edge-level posterior summaries
# -----------------------------------------------------------------------------

# Row-scaled factors Z (p x k x T) such that edge_ij(t) = sum_k Z[i,k,t] Z[j,k,t]
.edge_factors <- function(L, s2, type) {
  Tn <- dim(L)[3]; k <- dim(L)[2]
  for (t in seq_len(Tn)) {
    Lt <- matrix(L[, , t], ncol = k)
    if (type == "correlation") {
      # Sigma = L L' + s2 I
      L[, , t] <- Lt / sqrt(rowSums(Lt^2) + s2[t])
    } else {
      # Omega = (I - L A L') / s2, A = (s2 I + L'L)^-1  (Woodbury)
      # partial cor_ij = (L A L')_ij / sqrt((1 - (LAL')_ii)(1 - (LAL')_jj))
      R <- chol(solve(diag(s2[t], k) + crossprod(Lt)))
      Ft <- Lt %*% t(R)
      L[, , t] <- Ft / sqrt(1 - rowSums(Ft^2))
    }
  }
  L
}

#' Posterior summary of every edge
#'
#' @param bfit output of basil_fit().
#' @param type "correlation" (marginal; what BASIL's theory covers) or
#'   "partial" (partial correlations from the implied precision matrix;
#'   exploratory, see pipeline/README.md).
#' @param level credible level.
#' @param delta "practically non-zero" threshold for P(|edge| > delta | data).
#' @param features optional `features` data.frame to annotate assays.
#' @return list(edges = data.frame, mean/lower/upper = p x p matrices)
#'   Edge columns: mean, lower, upper, prob_direction = P(edge has the sign of
#'   its posterior mean), prob_gt_delta = P(|edge| > delta), bfdr = Bayesian FDR
#'   if all edges with prob_gt_delta at least this large are reported,
#'   ci_excludes_zero.
basil_edges <- function(bfit, type = c("correlation", "partial"), level = 0.95,
                        delta = 0.1, features = NULL, max_gb = 0.5) {
  .need("matrixStats")
  type <- match.arg(type)
  Z <- .edge_factors(bfit$Lambda_samples, bfit$sigma_sq_samples, type)
  p <- dim(Z)[1]; k <- dim(Z)[2]; Tn <- dim(Z)[3]
  a <- (1 - level) / 2
  M <- Lo <- Up <- Pd <- Pg <- matrix(NA_real_, p, p)
  b <- max(1, floor(max_gb * 1e9 / (8 * p * Tn)))
  for (i1 in seq(1, p - 1, by = b)) {
    rows <- i1:min(i1 + b - 1, p - 1); cols <- i1:p
    V <- array(0, c(length(rows), length(cols), Tn))
    for (t in seq_len(Tn))
      V[, , t] <- tcrossprod(matrix(Z[rows, , t], ncol = k), matrix(Z[cols, , t], ncol = k))
    dim(V) <- c(length(rows) * length(cols), Tn)
    q <- matrixStats::rowQuantiles(V, probs = c(a, 1 - a))
    put <- function(x) matrix(x, length(rows), length(cols))
    M[rows, cols] <- put(rowMeans(V))
    Lo[rows, cols] <- put(q[, 1]); Up[rows, cols] <- put(q[, 2])
    pp <- rowMeans(V > 0)
    Pd[rows, cols] <- put(pmax(pp, 1 - pp))
    Pg[rows, cols] <- put(rowMeans(abs(V) > delta))
  }
  ut <- which(upper.tri(M), arr.ind = TRUE)
  f <- bfit$features
  edges <- data.frame(feature1 = f[ut[, 1]], feature2 = f[ut[, 2]],
                      mean = M[ut], lower = Lo[ut], upper = Up[ut],
                      prob_direction = Pd[ut], prob_gt_delta = Pg[ut],
                      stringsAsFactors = FALSE)
  edges$ci_excludes_zero <- edges$lower > 0 | edges$upper < 0
  o <- order(-edges$prob_gt_delta, -abs(edges$mean))
  bf <- numeric(nrow(edges)); bf[o] <- cumsum(1 - edges$prob_gt_delta[o]) / seq_along(o)
  edges$bfdr <- bf
  if (!is.null(features)) {
    mod <- stats::setNames(features$modality, features$feature)
    edges$modality1 <- unname(mod[edges$feature1])
    edges$modality2 <- unname(mod[edges$feature2])
  }
  edges <- edges[o, ]; rownames(edges) <- NULL
  sym <- function(A) { A[lower.tri(A)] <- t(A)[lower.tri(A)]; dimnames(A) <- list(f, f); A }
  diag(M) <- diag(Lo) <- diag(Up) <- 1
  list(edges = edges, mean = sym(M), lower = sym(Lo), upper = sym(Up),
       type = type, level = level, delta = delta)
}

#' Pick the edges to report
#'
#' @param rule "ci": credible interval excludes zero (BASIL paper / vignette).
#'             "bfdr": P(|edge| > delta) large enough that the expected
#'             fraction of reported edges with |edge| <= delta is <= alpha.
#' @param min_abs additionally require |posterior mean| >= min_abs.
select_edges <- function(edges, rule = c("ci", "bfdr"), alpha = 0.05, min_abs = 0) {
  rule <- match.arg(rule)
  keep <- if (rule == "ci") edges$ci_excludes_zero else edges$bfdr <= alpha
  edges[keep & abs(edges$mean) >= min_abs, , drop = FALSE]
}

# -----------------------------------------------------------------------------
# 6. Plots
# -----------------------------------------------------------------------------

#' Forest plot: posterior mean and credible interval of the strongest edges
plot_edge_forest <- function(edges, top = 40, title = NULL) {
  .need("ggplot2")
  e <- edges[order(-abs(edges$mean)), ][seq_len(min(top, nrow(edges))), ]
  e$edge <- paste(e$feature1, "-", e$feature2)
  e$edge <- factor(e$edge, levels = rev(e$edge))
  e$sign <- ifelse(e$mean > 0, "positive", "negative")
  ggplot2::ggplot(e, ggplot2::aes(x = mean, y = edge, colour = sign)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_linerange(ggplot2::aes(xmin = lower, xmax = upper), linewidth = 0.6) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_colour_manual(values = c(positive = "#D73027", negative = "#4575B4")) +
    ggplot2::labs(x = "Edge strength (posterior mean, credible interval)", y = NULL,
                  title = title, colour = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Network of selected edges; nodes coloured by assay, edge width = |mean|,
#' edge opacity = how far the credible interval is from zero.
#' Only the `max_edges` strongest edges are drawn to keep the plot readable.
plot_basil_network <- function(edges, features = NULL, max_edges = 300, label_top = 40,
                               title = NULL, seed = 1) {
  .need("igraph"); .need("ggraph"); .need("ggplot2")
  if (nrow(edges) == 0) stop("No edges to plot.")
  e <- edges[order(-abs(edges$mean)), ][seq_len(min(max_edges, nrow(edges))), ]
  e$strength <- abs(e$mean)
  e$certainty <- pmin(abs(e$lower), abs(e$upper)) * !(e$lower < 0 & e$upper > 0)
  e$sign <- ifelse(e$mean > 0, "positive", "negative")
  g <- igraph::graph_from_data_frame(
    e[, c("feature1", "feature2", "strength", "certainty", "sign")], directed = FALSE)
  nm <- igraph::V(g)$name
  igraph::V(g)$assay <- if (is.null(features)) "feature" else
    features$modality[match(nm, features$feature)]
  st <- igraph::strength(g, weights = igraph::E(g)$strength)
  igraph::V(g)$node_strength <- st
  igraph::V(g)$label <- ifelse(rank(-st, ties.method = "first") <= label_top, nm, NA)
  set.seed(seed)
  ggraph::ggraph(g, layout = "fr", weights = igraph::E(g)$strength) +
    ggraph::geom_edge_link(ggplot2::aes(edge_width = strength, edge_alpha = certainty,
                                        edge_colour = sign)) +
    ggraph::scale_edge_width(range = c(0.2, 1.6), name = "|edge|") +
    ggraph::scale_edge_alpha(range = c(0.25, 0.95), name = "CI distance\nfrom 0") +
    ggraph::scale_edge_colour_manual(values = c(positive = "#D73027", negative = "#4575B4"),
                                     name = "Sign") +
    ggraph::geom_node_point(ggplot2::aes(colour = assay, size = node_strength)) +
    ggplot2::scale_size_continuous(range = c(1.5, 5), guide = "none") +
    ggplot2::scale_colour_brewer(palette = "Dark2", name = "Assay") +
    ggraph::geom_node_text(ggplot2::aes(label = label), repel = TRUE, size = 3,
                           na.rm = TRUE, max.overlaps = 30) +
    ggplot2::theme_void() + ggplot2::ggtitle(title)
}

#' Heatmap of posterior-mean edges, with non-selected edges set to zero
plot_edge_heatmap <- function(ed, selected, features = NULL, title = NULL) {
  .need("pheatmap")
  A <- ed$mean; A[] <- 0
  idx <- cbind(match(selected$feature1, rownames(A)), match(selected$feature2, rownames(A)))
  A[idx] <- selected$mean; A[idx[, 2:1, drop = FALSE]] <- selected$mean
  lim <- max(abs(A), 1e-8)
  ann <- if (is.null(features)) NA else
    data.frame(assay = features$modality, row.names = features$feature)
  pheatmap::pheatmap(A, breaks = seq(-lim, lim, length.out = 52),
                     color = grDevices::colorRampPalette(c("#4575B4", "white", "#D73027"))(51),
                     annotation_row = ann, annotation_col = ann, border_color = NA,
                     show_rownames = nrow(A) <= 60, show_colnames = nrow(A) <= 60,
                     treeheight_row = 0, treeheight_col = 0, main = title, silent = TRUE)
}

# -----------------------------------------------------------------------------
# One-call wrapper
# -----------------------------------------------------------------------------

#' Run the whole pipeline
#'
#' @param data,types,prep_args,covariates,inverse_normal see basil_prepare().
#' @param sets see basil_membership().
#' @param k,k_max,n_MC,seed,equalize_noise see basil_fit().
#' @param edge_type,level,delta see basil_edges().
#' @param rule,alpha,min_abs see select_edges().
#' @param out_dir if given, tables / plots / an .rds of the result are written there.
#' @return list(prep, C, fit, edges (all), selected, matrices, settings)
run_basil_pipeline <- function(data, types, sets = NULL, covariates = NULL,
                               prep_args = list(), inverse_normal = FALSE,
                               k = NULL, k_max = NULL, n_MC = 1000,
                               equalize_noise = TRUE,
                               edge_type = c("correlation", "partial"),
                               level = 0.95, delta = 0.1,
                               rule = c("ci", "bfdr"), alpha = 0.05, min_abs = 0.1,
                               out_dir = NULL, seed = 1, verbose = TRUE) {
  edge_type <- match.arg(edge_type); rule <- match.arg(rule)
  prep <- basil_prepare(data, types, prep_args, covariates, inverse_normal, verbose)
  C <- basil_membership(prep$features, sets, verbose = verbose)
  bfit <- basil_fit(prep$Y, C, k = k, k_max = k_max, n_MC = n_MC, seed = seed,
                    equalize_noise = equalize_noise, verbose = verbose)
  ed <- basil_edges(bfit, type = edge_type, level = level, delta = delta,
                    features = prep$features)
  sel <- select_edges(ed$edges, rule = rule, alpha = alpha, min_abs = min_abs)
  .msg(verbose, sprintf("Edges selected (%s): %d of %d", rule, nrow(sel), nrow(ed$edges)))
  if (length(unique(prep$features$modality)) > 1 && verbose && nrow(sel)) {
    tab <- table(paste(pmin(sel$modality1, sel$modality2),
                       pmax(sel$modality1, sel$modality2), sep = " - "))
    for (i in seq_along(tab)) message(sprintf("  %-28s %d", names(tab)[i], tab[i]))
  }
  res <- list(prep = prep, C = C, fit = bfit, edges = ed$edges, selected = sel,
              matrices = ed[c("mean", "lower", "upper")],
              settings = list(edge_type = edge_type, level = level, delta = delta,
                              rule = rule, alpha = alpha, min_abs = min_abs,
                              k = bfit$k, n_MC = n_MC, seed = seed))
  if (!is.null(out_dir)) {
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(ed$edges, file.path(out_dir, "edges_all.csv"), row.names = FALSE)
    utils::write.csv(sel, file.path(out_dir, "edges_selected.csv"), row.names = FALSE)
    utils::write.csv(prep$prep_log, file.path(out_dir, "feature_prep_log.csv"), row.names = FALSE)
    slim <- res; slim$fit$Lambda_samples <- NULL   # draws can be large; refit with same seed
    saveRDS(slim, file.path(out_dir, "basil_result.rds"))
    lab <- sprintf("BASIL %s network (%g%% CI)", edge_type, 100 * level)
    if (nrow(sel)) {
      ggplot2::ggsave(file.path(out_dir, "edge_forest.pdf"),
                      plot_edge_forest(sel, title = lab), width = 8,
                      height = 2 + 0.2 * min(40, nrow(sel)))
      ggplot2::ggsave(file.path(out_dir, "network.pdf"),
                      plot_basil_network(sel, prep$features, title = paste0(
                        lab, if (nrow(sel) > 300) ", 300 strongest edges" else "")),
                      width = 10, height = 8)
      grDevices::pdf(file.path(out_dir, "edge_heatmap.pdf"), width = 9, height = 8)
      grid::grid.newpage()
      grid::grid.draw(plot_edge_heatmap(ed, sel, prep$features, title = lab)$gtable)
      grDevices::dev.off()
    }
    .msg(verbose, "Results written to ", normalizePath(out_dir))
  }
  invisible(res)
}
