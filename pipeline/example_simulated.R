# Worked example on simulated data shaped like the real assays.
# Run from the repository root:  Rscript pipeline/example_simulated.R
source("pipeline/basil_pipeline.R")
set.seed(2026)

# ---- simulate 120 subjects, 4 shared latent factors -------------------------
n <- 120; k <- 4
p <- c(bile_acids = 24, scfa = 8, olink = 92, rna = 150)
eta <- matrix(rnorm(n * k), n, k)
latent <- function(pj) {                       # sparse loadings + noise
  L <- matrix(rnorm(pj * k) * rbinom(pj * k, 1, 0.4), pj, k)
  eta %*% t(L) + matrix(rnorm(n * pj), n, pj)
}
ids <- sprintf("S%03d", 1:n)

ba <- 2^(latent(p["bile_acids"]) + 3)                       # concentrations
colnames(ba) <- paste0("BA", 1:ncol(ba))
ba[ba < quantile(ba, 0.08)] <- NA                           # below LOQ -> NA
scfa <- 2^(latent(p["scfa"]) + 6); colnames(scfa) <- paste0("SCFA", 1:ncol(scfa))
npx <- latent(p["olink"]) + 5; colnames(npx) <- paste0("PROT", 1:ncol(npx))
mu <- 2^(latent(p["rna"]) * 0.6 + 6)
rna <- matrix(rnbinom(length(mu), mu = mu, size = 20), n); colnames(rna) <- paste0("GENE", 1:ncol(rna))
rownames(ba) <- rownames(scfa) <- rownames(npx) <- rownames(rna) <- ids
meta <- data.frame(age = rnorm(n, 55, 10), sex = factor(sample(c("F", "M"), n, TRUE)),
                   row.names = ids)

# ---- optional prior knowledge: feature sets (may span assays) ---------------
sets <- list(
  primary_BA      = paste0("BA", 1:8),
  secondary_BA    = paste0("BA", 9:16),
  conjugated_BA   = paste0("BA", c(1:4, 9:12, 17:20)),
  inflammation    = c(paste0("PROT", 1:20), paste0("GENE", 1:25)),
  barrier         = c(paste0("PROT", 21:35), paste0("GENE", 26:60))
)

# ---- run --------------------------------------------------------------------
res <- run_basil_pipeline(
  data  = list(bile_acids = ba, scfa = scfa, olink = npx, rna = rna),
  types = c(bile_acids = "concentration", scfa = "concentration",
            olink = "npx", rna = "counts"),
  sets = sets,
  covariates = meta,
  out_dir = "pipeline/example_output"
)

head(res$selected, 10)
# cross-assay edges only, strongest first
cross <- subset(res$selected, modality1 != modality2)
head(cross[order(-abs(cross$mean)), ], 10)
