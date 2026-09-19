# Write an ADMIXTURE warm-start Q from the fastICA structure solution, so that
# ADMIXTURE can be started at the x|x| answer instead of at random.
#
# ADMIXTURE reads initial values from <prefix>.<K>.Q.in (and optionally
# <prefix>.<K>.P.in) sitting next to the .bed. We supply Q only and let it
# estimate P: our .bed codes the derived allele as A2 while the .bim allele
# labels come from REF/ALT, so which allele ADMIXTURE's P refers to is not
# something we should guess at. Q is unaffected by that.
#
# The memberships are the same ones the structure plots use: A = pmax(L, 0)
# for the K = 40 rank-1 maxima, with columns rescaled so the rows of B sum to
# 1. Those are per haplotype, so each row is normalized to sum to 1 and the
# two haplotypes of an individual are averaged, giving a Q whose rows sum to 1.
#
# Usage: Rscript code/write_admixture_init.R <outdir> <prefix>
# Writes <outdir>/<prefix>.<K>.Q.in

args   <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]; prefix <- args[2]

meta   <- readRDS("data/hgdp_meta.rds")
fit    <- readRDS("output/hgdp_fastica.rds")
struct <- readRDS("output/hgdp_structure.rds")

L <- fit$fits$k40$maxima$L
r <- struct$maxima_k40$r
hm <- meta[match(sub("_[12]$", "", fit$hap_ids), meta$sample), ]

# same ordering and scaling as the structure plots
w <- pmax(L, 0)^2
o <- order(colSums(w * as.integer(hm$pop)) / colSums(w))
A <- sweep(pmax(L[, o, drop = FALSE], 0), 2, r[o], "*")

stopifnot(all(rowSums(A) > 0))
P <- A / rowSums(A)                      # per-haplotype memberships
Q <- 0.5 * (P[c(TRUE, FALSE), , drop = FALSE] +
            P[c(FALSE, TRUE), , drop = FALSE])   # average the two chromosomes

# ADMIXTURE does not like exact zeros in the initial Q; floor and renormalize.
eps <- 1e-5
Q <- pmax(Q, eps)
Q <- Q / rowSums(Q)

ids <- unique(sub("_[12]$", "", fit$hap_ids))
stopifnot(nrow(Q) == length(ids), identical(ids, as.character(meta$sample)))
K <- ncol(Q)

f <- file.path(outdir, sprintf("%s.%d.Q.in", prefix, K))
write.table(format(Q, digits = 6, scientific = FALSE), f,
            quote = FALSE, row.names = FALSE, col.names = FALSE)
cat(sprintf("wrote %s: %d individuals x %d components\n", f, nrow(Q), K))
cat("row sums range:", range(rowSums(Q)), "\n")
cat("column means (geographic order):\n"); print(round(colMeans(Q), 3))
