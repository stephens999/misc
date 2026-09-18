# Is there structure left in the residuals of the non-negative fit?
#
# Takes the K = 40 maxima fit X ~ A B of code/fit_hgdp_structure.R, with
# A = pmax(L, 0) and B the least-squares coefficients, forms the residual
# E = X - A B, and runs the same rank-1 fastICA search on it that the main
# analysis ran on X: centre, whiten to 40 dimensions, 1000 independent random
# starts, cluster the results. This is the HGDP counterpart of
# code/fit_tgp_residual.R.
#
# E is 1858 x 215468 and never needs to be formed. Its Gram matrix is
#   E E' = G_X - M A' - A M' + A (B B') A',   M = X B',
# and the centred version follows from E E' and mu_E = colMeans(E) as usual.
# G_X, M, A'X and the column means are all accumulated over SNP blocks, so no
# numeric copy of the genotype matrix is held.
#
# It also saves a few per-haplotype and per-population summaries that the
# analysis uses to interpret the residual (centred norms, the residual share
# of each, population allele-frequency deviation from the global mean, and
# heterozygosity), since they need the same pass over the genotypes.
#
# Usage: Rscript code/fit_hgdp_residual.R
# Writes output/hgdp_residual.rds.

K_res  <- 40     # whitening dimension for the residual search
n_r1   <- 1000
n_iter <- 500

f   <- readRDS("output/hgdp_fastica.rds")
dat <- readRDS("data/hgdp_hap_matrix.rds")
Y   <- dat$Y
n   <- nrow(Y)
p   <- ncol(Y)

blk    <- 20000
blocks <- split(seq_len(p), ceiling(seq_len(p) / blk))
as_num <- function(cols) matrix(as.numeric(Y[, cols, drop = FALSE]) - 48, nrow = n)

A  <- pmax(f$fits$k40$maxima$L, 0)
stopifnot(nrow(A) == n)
cat(sprintf("%d haplotypes x %d SNPs, %d sources\n", n, p, ncol(A)))

# ---- B, the Gram matrix of X, and X B', in one pass -------------------------

AtX <- matrix(0, ncol(A), p)
G_X <- matrix(0, n, n)
mu_X <- numeric(p)
for (cols in blocks) {
  Yb          <- as_num(cols)
  AtX[, cols] <- crossprod(A, Yb)
  G_X         <- G_X + tcrossprod(Yb)
  mu_X[cols]  <- colMeans(Yb)
}
G_X <- (G_X + t(G_X)) / 2
B   <- solve(crossprod(A), AtX)
cat("fraction of B below -0.01:", signif(mean(B < -0.01), 3),
    " min B:", round(min(B), 4), "\n")

# M = X B', accumulated the same way. The same loop collects the per-haplotype
# centred norms, the per-haplotype centred residual norms, and the population
# allele frequencies, all of which the analysis needs.
meta   <- readRDS("data/hgdp_meta.rds")
hpop   <- meta$pop[match(sub("_[12]$", "", f$hap_ids), meta$sample)]
pops   <- levels(hpop)
mu_E0  <- mu_X - as.vector(colMeans(A) %*% B)

M   <- matrix(0, n, ncol(A))
cx2 <- numeric(n)                       # || centred X row ||^2
ce2 <- numeric(n)                       # || centred residual row ||^2
Fq  <- matrix(0, length(pops), p, dimnames = list(pops, NULL))
for (cols in blocks) {
  Yb  <- as_num(cols)
  M   <- M + Yb %*% t(B[, cols, drop = FALSE])
  cx2 <- cx2 + rowSums(sweep(Yb, 2, mu_X[cols], "-")^2)
  Eb  <- sweep(Yb - A %*% B[, cols, drop = FALSE], 2, mu_E0[cols], "-")
  ce2 <- ce2 + rowSums(Eb^2)
  for (k in seq_along(pops))
    Fq[k, cols] <- colMeans(Yb[hpop == pops[k], , drop = FALSE])
}
pop_dev <- rowSums(sweep(Fq, 2, mu_X, "-")^2)   # || f_pop - mu ||^2
pop_het <- rowMeans(2 * Fq * (1 - Fq))
rm(Fq); gc()

# ---- Gram matrix of the residual, and its centred version -------------------

G_E <- G_X - tcrossprod(M, A) - tcrossprod(A, M) + A %*% tcrossprod(B) %*% t(A)
G_E <- (G_E + t(G_E)) / 2

ss_X <- sum(diag(G_X)); ss_E <- sum(diag(G_E))
cat("residual sum of squares / total:", round(ss_E / ss_X, 4), "\n")

# mu_E = colMeans(X) - colMeans(A) B, and a_E = E mu_E = X mu_E - A (B mu_E)
mu_E <- mu_X - as.vector(colMeans(A) %*% B)
a_E  <- numeric(n)
for (cols in blocks) a_E <- a_E + as_num(cols) %*% mu_E[cols]
a_E  <- as.vector(a_E) - as.vector(A %*% (B %*% mu_E))

G_C <- G_E - outer(a_E, rep(1, n)) - outer(rep(1, n), a_E) + sum(mu_E^2)
G_C <- (G_C + t(G_C)) / 2
rm(Y, G_X, AtX); gc()

eig <- eigen(G_C, symmetric = TRUE)
d   <- sqrt(pmax(eig$values, 0))
U_c <- eig$vectors[, 1:K_res, drop = FALSE]
rownames(U_c) <- f$hap_ids
cat("top singular values of the centred residual:", round(d[1:12]), "\n")

# ---- rank-1 fastICA on the residual -----------------------------------------

r1_update <- function(U, W) {
  P <- t(U) %*% W
  W <- U %*% (2 * abs(P)) - sweep(W, 2, colSums(2 * sign(P)), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

U <- sqrt(n) * t(U_c)
set.seed(1)
W <- matrix(rnorm(K_res * n_r1), K_res, n_r1)
W <- sweep(W, 2, sqrt(colSums(W^2)), "/")
for (i in seq_len(n_iter)) { W_old <- W; W <- r1_update(U, W) }
converged <- max(1 - abs(colSums(W * W_old)))

L  <- t(U) %*% W
cl <- cutree(hclust(as.dist(1 - cor(L)), method = "complete"), h = 0.05)
j  <- match(unique(cl), cl)
L_max <- L[, j, drop = FALSE]
rownames(L_max) <- f$hap_ids
size  <- as.integer(table(cl)[as.character(cl[j])])
cat("distinct maxima in the residual:", length(j), "| basins:",
    paste(sort(size, decreasing = TRUE), collapse = ","), "\n")

saveRDS(list(K = K_res, d = d[1:60], U_c = U_c, ids = f$hap_ids,
             cx2 = cx2, ce2 = ce2, pop_dev = pop_dev, pop_het = pop_het,
             L = L_max, size = size, converged = converged,
             ss_ratio = ss_E / ss_X, n_sources = ncol(A),
             n_cut = sapply(c(0.90, 0.95, 0.99), function(tau)
               length(unique(cutree(hclust(as.dist(1 - cor(L)),
                                           method = "complete"), h = 1 - tau))))),
        "output/hgdp_residual.rds")
cat("wrote output/hgdp_residual.rds\n")
