# Scaling needed to turn the fastICA sources of code/fit_tgp_fastica.R into
# non-negative memberships for structure plots.
#
# For a source matrix L we set A = pmax(L, 0) -- the x|x| sources are one-sided,
# so their positive parts carry the signal -- and fit
#
#   X ~ A B + E
#
# by least squares, with no intercept column. If B is non-negative we can scale
# each row of B to sum to 1 and scale the corresponding column of A to
# compensate, leaving the rows of A as memberships on a common scale. Whether
# B >= 0 actually holds is an empirical question and is the main thing this
# script reports.
#
# Only rowSums(B) is needed to do the rescaling, so that is what we save rather
# than B itself, which would be K x 185116 per fit.
#
# Usage: Rscript code/fit_tgp_structure.R
# Reads output/tgp_fastica.rds and data/tgp_data_matrix.rds.
# Writes output/tgp_structure.rds.

f <- readRDS("output/tgp_fastica.rds")
Y <- readRDS("data/tgp_data_matrix.rds")
n <- nrow(Y)

chrom <- as.integer(sub(":.*", "", colnames(Y)))
snp_sets <- list(all  = seq_len(ncol(Y)),
                 odd  = which(chrom %% 2 == 1),
                 even = which(chrom %% 2 == 0))

# Least-squares B for given sources, over a subset of individuals and SNPs.
# `keep` indexes rows of Y, `cols` indexes columns.
ab_fit <- function(L, keep, cols) {
  A   <- pmax(L, 0)
  Ap  <- matrix(0, n, ncol(A))
  Ap[keep, ] <- A
  AtX <- crossprod(Ap, Y[, cols, drop = FALSE])
  B   <- solve(crossprod(A), AtX)
  list(K        = ncol(A),
       n_ind    = length(keep),
       n_snps   = length(cols),
       r        = rowSums(pmax(B, 0)),   # rescales columns of A
       # Correlations between the components' estimated allele frequencies,
       # across SNPs. Invariant to the per-row rescaling above, so this is the
       # same whether computed from B or from the row-normalized version.
       cor_B    = cor(t(B)),
       # Also the K x K cross-product and row means, which allow the same
       # comparison after removing the frequency spectrum that all components
       # share: P (BB') P with P = I - 11'/K.
       gram_B   = tcrossprod(B),
       mean_B   = rowMeans(B),
       frac_neg = mean(B < 0),
       min_B    = min(B),
       frac_lt  = mean(B < -0.01),
       mean_abs = mean(abs(B)),
       r2       = sum(B * AtX) / sum(Y[keep, cols, drop = FALSE]^2))
}

# Stability pruning of the K = 30 rank-r fits: match the odd-half factors to the
# even-half ones so as to maximize total (signed) correlation, and keep only the
# pairs that agree. The matched correlations split cleanly into 15 above 0.48
# and 15 below 0.09, so the cutoff is not delicate.
C_half  <- cor(f$split$odd$fit$L, f$split$even$fit$L)
perm    <- local({
  pr <- RcppHungarian::HungarianSolver(-C_half)$pairs
  pr[order(pr[, 1]), 2]
})
match_cor <- C_half[cbind(seq_along(perm), perm)]
stable    <- which(match_cor > 0.3)
cat("replicating K = 30 factors:", length(stable), "of", ncol(C_half), "\n")

# Every set of sources saved by fit_tgp_fastica.R, each against the SNPs it was
# estimated from: the chromosome halves are fit on their own half.
spec <- list(
  full_k20        = list(L = f$fits$k20$L,          ids = f$ids$full,      snps = "all"),
  full_k30        = list(L = f$fits$k30$L,          ids = f$ids$full,      snps = "all"),
  full_maxima     = list(L = f$maxima$full$L,       ids = f$ids$full,      snps = "all"),
  pruned_k30      = list(L = f$no_sparse$fit$L,     ids = f$ids$no_sparse, snps = "all"),
  pruned_maxima   = list(L = f$maxima$no_sparse$L,  ids = f$ids$no_sparse, snps = "all"),
  odd_k30         = list(L = f$split$odd$fit$L,     ids = f$ids$odd,       snps = "odd"),
  odd_maxima      = list(L = f$maxima$odd$L,        ids = f$ids$odd,       snps = "odd"),
  even_k30        = list(L = f$split$even$fit$L,    ids = f$ids$even,      snps = "even"),
  even_maxima     = list(L = f$maxima$even$L,       ids = f$ids$even,      snps = "even"),
  odd_k30_stable  = list(L = f$split$odd$fit$L[, stable, drop = FALSE],
                         ids = f$ids$odd,  snps = "odd"),
  even_k30_stable = list(L = f$split$even$fit$L[, perm[stable], drop = FALSE],
                         ids = f$ids$even, snps = "even"))

# The K = 14 B itself is saved separately, for the pairwise scatter plots of
# component allele frequencies. 14 x 185116 is too big to keep for every fit,
# but worth it for the one the analysis focuses on.
save_B_for <- "pruned_maxima"

structure_fits <- lapply(names(spec), function(nm) {
  sp   <- spec[[nm]]
  keep <- match(sp$ids, rownames(Y))
  stopifnot(!any(is.na(keep)), nrow(sp$L) == length(keep))
  out  <- ab_fit(sp$L, keep, snp_sets[[sp$snps]])
  if (identical(nm, save_B_for)) {
    A   <- pmax(sp$L, 0)
    Ap  <- matrix(0, n, ncol(A)); Ap[keep, ] <- A
    cols <- snp_sets[[sp$snps]]
    B   <- solve(crossprod(A), crossprod(Ap, Y[, cols, drop = FALSE]))
    colnames(B) <- colnames(Y)[cols]
    saveRDS(B, sprintf("output/tgp_B_%s.rds", nm))
  }
  cat(sprintf("%-15s K=%2d  frac B<0 %.4f  min B %+.4f  frac<-0.01 %.4f  R2 %.3f\n",
              nm, out$K, out$frac_neg, out$min_B, out$frac_lt, out$r2))
  c(out, list(snps = sp$snps))
})
names(structure_fits) <- names(spec)

# The matching is saved too, so the analysis file can label the stable factors
# and report which of the 30 were kept.
saveRDS(c(structure_fits,
          list(half_match = list(perm = perm, match_cor = match_cor,
                                 stable = stable))),
        "output/tgp_structure.rds")
