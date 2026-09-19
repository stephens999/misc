# Is there structure left in the residuals of the non-negative fit?
#
# Takes the K = 14 stable-maxima fit X ~ A B from code/fit_tgp_structure.R,
# forms the residual E = X - A B, and runs the same rank-1 fastICA search on it
# that the main analysis ran on X: centre, whiten to 30 dimensions, 1000
# independent random starts, cluster the results.
#
# E is 2464 x 185116 and never needs to be formed. Its Gram matrix is
#   E E' = G_X - M A' - A M' + A (B B') A',    M = X B',
# and the centred version follows from E E' and mu_E = colMeans(E) as usual.
#
# Usage: Rscript code/fit_tgp_residual.R
# Writes output/tgp_residual.rds.

f <- readRDS("output/tgp_fastica.rds")
B <- readRDS("output/tgp_B_pruned_maxima.rds")
Y <- readRDS("data/tgp_data_matrix.rds")

ids  <- f$ids$no_sparse
keep <- match(ids, rownames(Y))
A    <- pmax(f$maxima$no_sparse$L, 0)       # columns match the rows of B
nk   <- nrow(A)
stopifnot(nrow(A) == length(keep), ncol(A) == nrow(B))

# ---- Gram matrix of the residual, without forming it ------------------------

Ap  <- matrix(0, nrow(Y), ncol(A)); Ap[keep, ] <- A
M   <- crossprod(t(Y[keep, , drop = FALSE]), t(B))      # X B' : nk x K
G_X <- tcrossprod(Y[keep, , drop = FALSE])
G_E <- G_X - tcrossprod(M, A) - tcrossprod(A, M) + A %*% tcrossprod(B) %*% t(A)
G_E <- (G_E + t(G_E)) / 2

ss_X <- sum(diag(G_X)); ss_E <- sum(diag(G_E))
cat("residual sum of squares / total:", round(ss_E / ss_X, 4), "\n")

# centre: mu_E = colMeans(X) - colMeans(A) B
mu_E <- colMeans(Y[keep, , drop = FALSE]) - as.vector(colMeans(A) %*% B)
a_E  <- as.vector(Y[keep, , drop = FALSE] %*% mu_E) - as.vector(A %*% (B %*% mu_E))
G_C  <- G_E - outer(a_E, rep(1, nk)) - outer(rep(1, nk), a_E) + sum(mu_E^2)
G_C  <- (G_C + t(G_C)) / 2
rm(Y, G_X, Ap); gc()

eig <- eigen(G_C, symmetric = TRUE)
d   <- sqrt(pmax(eig$values, 0))
U_c <- eig$vectors[, 1:30, drop = FALSE]
rownames(U_c) <- ids
cat("top singular values of the centred residual:", round(d[1:12]), "\n")

# ---- rank-1 fastICA on the residual -----------------------------------------

r1_update <- function(U, W) {
  P <- t(U) %*% W
  W <- U %*% (2 * abs(P)) - sweep(W, 2, colSums(2 * sign(P)), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

U <- sqrt(nk) * t(U_c)
set.seed(1)
W <- matrix(rnorm(30 * 1000), 30, 1000)
W <- sweep(W, 2, sqrt(colSums(W^2)), "/")
for (i in 1:500) W <- r1_update(U, W)

L  <- t(U) %*% W
cl <- cutree(hclust(as.dist(1 - cor(L)), method = "complete"), h = 0.05)
j  <- match(unique(cl), cl)
L_max <- L[, j, drop = FALSE]
rownames(L_max) <- ids
size  <- as.integer(table(cl)[as.character(cl[j])])
cat("distinct maxima in the residual:", length(j), "| basins:",
    paste(sort(size, decreasing = TRUE), collapse = ","), "\n")

saveRDS(list(d = d, U_c = U_c, ids = ids, L = L_max, size = size,
             ss_ratio = ss_E / ss_X),
        "output/tgp_residual.rds")
