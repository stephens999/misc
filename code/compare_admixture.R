# Compare the two ADMIXTURE runs at K = 25 with each other and with the
# fastICA memberships they are meant to be compared against:
#
#   * output/admixture/hgdp/       random start (--seed=1)
#   * output/admixture/hgdp_init/  started at the normalized fastICA solution
#     (-Q, reading hgdp.25.Q.in written by code/write_admixture_init.R)
#
# Both use the same 929 x 215468 .bed, so their log-likelihoods are directly
# comparable: a higher one means the warm start reached a better optimum.
#
# Components are matched across solutions by maximizing total correlation
# (Hungarian assignment) so the three structure plots can share a colour.
#
# Usage: Rscript code/compare_admixture.R
# Writes output/hgdp_admixture_compare.rds.

suppressMessages({library(RcppHungarian)})

meta   <- readRDS("data/hgdp_meta.rds")
fit    <- readRDS("output/hgdp_fastica.rds")
struct <- readRDS("output/hgdp_structure.rds")

dir_rand <- "output/admixture/hgdp_thin"
dir_init <- "output/admixture/hgdp_thin_init"

# Log-likelihood, block-iteration count and compute time from an ADMIXTURE log.
run_info <- function(f) {
  x  <- readLines(f)
  ll <- grep("^Loglikelihood:", x, value = TRUE)
  it <- grep("QN/Block", x, value = TRUE)
  secs <- as.numeric(sub(".*Elapsed: ([0-9.]+).*", "\\1", it))
  c(loglik = if (!length(ll)) NA_real_ else
      as.numeric(sub("^Loglikelihood:\\s*", "", tail(ll, 1))),
    iterations = length(it), minutes = sum(secs) / 60)
}
info <- rbind(random  = run_info(file.path(dir_rand, "admixture_k25.log")),
              fastICA = run_info(file.path(dir_init, "admixture_k25_init.log")))
ll <- info[, "loglik"]

# Any extra random seeds that have finished (output/admixture/seeds/s*/).
seed_dirs <- Sys.glob("output/admixture/seeds/s*/admixture.log")
seeds <- if (length(seed_dirs)) {
  z <- t(sapply(seed_dirs, run_info))
  rownames(z) <- sub(".*/s([0-9]+)/.*", "seed \\1", seed_dirs)
  z[!is.na(z[, "loglik"]), , drop = FALSE]
} else NULL

Qr <- as.matrix(read.table(file.path(dir_rand, "hgdp.25.Q")))
Qi <- as.matrix(read.table(file.path(dir_init, "hgdp.25.Q")))
ids <- read.table(file.path(dir_rand, "hgdp.fam"))$V1
stopifnot(nrow(Qr) == length(ids), identical(as.character(ids),
                                             as.character(meta$sample)))

# fastICA memberships, per individual, in the same form (rows sum to 1)
L  <- fit$fits$k40$maxima$L
hm <- meta[match(sub("_[12]$", "", fit$hap_ids), meta$sample), ]
w  <- pmax(L, 0)^2
o  <- order(colSums(w * as.integer(hm$pop)) / colSums(w))
A  <- sweep(pmax(L[, o, drop = FALSE], 0), 2, struct$maxima_k40$r[o], "*")
Pm <- A / rowSums(A)
Qf <- 0.5 * (Pm[c(TRUE, FALSE), , drop = FALSE] + Pm[c(FALSE, TRUE), , drop = FALSE])

# Match the two ADMIXTURE solutions onto the fastICA column order.
match_to <- function(X, ref) {
  C <- cor(ref, X)
  pr <- HungarianSolver(-C)$pairs
  perm <- pr[order(pr[, 1]), 2]
  list(X = X[, perm, drop = FALSE], r = C[cbind(seq_along(perm), perm)])
}
mr <- match_to(Qr, Qf); mi <- match_to(Qi, Qf)
mri <- match_to(mi$X, mr$X)      # the two ADMIXTURE runs against each other

pop_label <- function(q, m, top = 2) {
  v <- sort(tapply(q, m$pop, mean), decreasing = TRUE)
  paste(names(v)[seq_len(top)], collapse = "/")
}
ipop <- meta[, c("sample", "pop", "region")]

out <- list(
  loglik  = ll,
  info    = info,
  seeds   = seeds,
  ids     = as.character(ids),
  pop     = meta$pop,
  region  = meta$region,
  Q_fastica = Qf,
  Q_random  = mr$X, r_random  = mr$r,
  Q_init    = mi$X, r_init    = mi$r,
  r_between = mri$r,
  labels    = apply(Qf, 2, pop_label, m = ipop))
saveRDS(out, "output/hgdp_admixture_compare.rds")

cat("=== log-likelihoods (same data, same K; higher is better) ===\n")
print(round(info, 2))
if (!is.null(seeds)) { cat("\nadditional random seeds:\n"); print(round(seeds, 2)) }
cat("difference (fastICA-init - random):", round(ll[2] - ll[1], 2), "\n\n")
cat("=== component matching, correlation with the fastICA memberships ===\n")
print(data.frame(component = out$labels, random = round(mr$r, 2),
                 fastICA_init = round(mi$r, 2), row.names = NULL))
cat("\nmedian r to fastICA: random", round(median(mr$r), 2),
    " init", round(median(mi$r), 2), "\n")
cat("median r between the two ADMIXTURE runs:", round(median(mri$r), 2), "\n")
cat("components matching between the runs at r > 0.99:", sum(mri$r > 0.99),
    "of", ncol(Qr), "\n")
