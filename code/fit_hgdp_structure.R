# Scaling needed to turn the fastICA sources of code/fit_hgdp_fastica.R into
# non-negative memberships for structure plots. This is the HGDP counterpart
# of code/fit_tgp_structure.R and does the same thing.
#
# For a source matrix L we set A = pmax(L, 0) -- the x|x| sources are one-sided,
# so their positive parts carry the signal -- and fit
#
#   X ~ A B + E
#
# by least squares, with no intercept column and on the uncentered X. If B is
# non-negative we can scale each row of B to sum to 1 and scale the
# corresponding column of A to compensate, leaving the rows of A as
# memberships on a common scale. Whether B >= 0 actually holds is an empirical
# question and is the main thing this script reports.
#
# Only rowSums(pmax(B, 0)) is needed for that rescaling, so that is what we
# save rather than B itself, which would be K x 215468 per fit.
#
# X here is the haplotype matrix, so the rows of A are haplotypes rather than
# individuals and a "membership" is a membership of one chromosome.
#
# Usage: Rscript code/fit_hgdp_structure.R
# Reads output/hgdp_fastica.rds and data/hgdp_hap_matrix.rds.
# Writes output/hgdp_structure.rds.

f   <- readRDS("output/hgdp_fastica.rds")
dat <- readRDS("data/hgdp_hap_matrix.rds")
Y   <- dat$Y                      # raw, nhap x nsnp, entries 0x30 / 0x31
n   <- nrow(Y)
p   <- ncol(Y)

blk    <- 20000
blocks <- split(seq_len(p), ceiling(seq_len(p) / blk))
as_num <- function(cols) matrix(as.numeric(Y[, cols, drop = FALSE]) - 48, nrow = n)

# Least-squares B for given sources over a subset of the haplotypes. A'X and
# the total sum of squares are accumulated over SNP blocks so that no numeric
# copy of the full matrix is ever formed.
ab_fit <- function(L, keep, save_B = FALSE) {
  A  <- pmax(L, 0)
  Ap <- matrix(0, n, ncol(A))
  Ap[keep, ] <- A

  AtX <- matrix(0, ncol(A), p)
  ssq <- 0
  for (cols in blocks) {
    Yb           <- as_num(cols)
    AtX[, cols]  <- crossprod(Ap, Yb)
    ssq          <- ssq + sum(Yb[keep, , drop = FALSE]^2)
  }
  B <- solve(crossprod(A), AtX)
  if (save_B) colnames(B) <- colnames(Y)

  out <- list(K        = ncol(A),
              n_hap    = length(keep),
              n_snps   = p,
              r        = rowSums(pmax(B, 0)),   # rescales the columns of A
              cor_B    = cor(t(B)),
              frac_neg = mean(B < 0),
              min_B    = min(B),
              frac_lt  = mean(B < -0.01),
              n_lt     = sum(B < -0.01),
              mean_abs = mean(abs(B)),
              r2       = sum(B * AtX) / ssq)
  if (save_B) out$B <- B
  out
}

# The source sets worth comparing: at each whitening dimension, the rank-1
# maxima and the rank-r factors.
spec <- list(
  maxima_k40 = list(L = f$fits$k40$maxima$L, ids = f$hap_ids),
  rankr_k40  = list(L = f$fits$k40$fit$L,    ids = f$hap_ids),
  maxima_k60 = list(L = f$fits$k60$maxima$L, ids = f$hap_ids),
  rankr_k60  = list(L = f$fits$k60$fit$L,    ids = f$hap_ids))

structure_fits <- lapply(names(spec), function(nm) {
  sp   <- spec[[nm]]
  keep <- match(sp$ids, rownames(Y))
  stopifnot(!any(is.na(keep)), nrow(sp$L) == length(keep))
  out  <- ab_fit(sp$L, keep)
  cat(sprintf("%-14s K=%2d  frac B<0 %.4f  min B %+.4f  n<-0.01 %6d of %9d  mean|B| %.3f  R2 %.3f\n",
              nm, out$K, out$frac_neg, out$min_B, out$n_lt,
              out$K * out$n_snps, out$mean_abs, out$r2))
  out
})
names(structure_fits) <- names(spec)

saveRDS(structure_fits, "output/hgdp_structure.rds")
cat("wrote output/hgdp_structure.rds\n")
