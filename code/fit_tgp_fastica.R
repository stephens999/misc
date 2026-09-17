# Run rank-r fastICA with the x|x| contrast and polar orthogonalization on the
# 1000 Genomes Project phase 3 genotype matrix. The driver is the one from
# single-cell-jamboree/analysis/pancreas_celseq2_ica_02.Rmd.
#
# The genotype matrix (2504 individuals x 185116 SNPs, counts of the derived
# allele) is ~3.7 Gb as doubles, so we never form a centred copy of it. Instead
# we get the left singular vectors from the centred Gram matrix
#
#   C C' = Y Y' - a 1' - 1 a' + (mu'mu) 11',   a = Y mu,  mu = colMeans(Y),
#
# which is only 2504 x 2504. The whitened row space U_w is then the same input
# that the jamboree analysis feeds to fastICA.
#
# Usage: Rscript code/fit_tgp_fastica.R
# Writes output/tgp_fastica.rds.

K        <- c(20, 30)    # whitening dimensions to fit (r = n.comp in each)
n_iter   <- 200
n_starts <- 5            # random starts per K; we keep the best objective

Y <- readRDS("data/tgp_data_matrix.rds")
n <- nrow(Y)

# ---- centred SVD via the Gram matrix -----------------------------------------

G_raw  <- tcrossprod(Y)
csums  <- colSums(Y)

# Centred Gram matrix for a subset of the rows, built from the raw Gram matrix
# so that we never form a second copy of Y. With mu the column means of the
# subset and a = Y[keep, ] mu,
#   C C' = (YY')[keep, keep] - a 1' - 1 a' + (mu'mu) 11'.
centred_gram <- function(keep) {
  nk <- length(keep)
  mu <- (csums - colSums(Y[setdiff(seq_len(n), keep), , drop = FALSE])) / nk
  a  <- as.vector(Y[keep, , drop = FALSE] %*% mu)
  G  <- G_raw[keep, keep, drop = FALSE] -
          outer(a, rep(1, nk)) - outer(rep(1, nk), a) + sum(mu^2)
  (G + t(G)) / 2   # symmetrize away rounding asymmetry
}

eigen_of <- function(keep) {
  e <- eigen(centred_gram(keep), symmetric = TRUE)
  U <- e$vectors
  rownames(U) <- rownames(Y)[keep]
  list(d = sqrt(pmax(e$values, 0)), U_c = U)   # d: singular values of centred Y
}

full   <- eigen_of(seq_len(n))
d      <- full$d
U_c    <- full$U_c

# Note: Y and G_raw are kept alive until after the fits below, because the
# second, reduced-data SVD depends on which individuals the fits single out.

# ---- fastICA -----------------------------------------------------------------

# Polar factor via eigendecomposition (avoids SVD column reordering).
polar <- function(W) {
  eig      <- eigen(t(W) %*% W, symmetric = TRUE)
  Ainvhalf <- eig$vectors %*%
               diag(1 / sqrt(pmax(eig$values, 1e-14))) %*%
               t(eig$vectors)
  W %*% Ainvhalf
}

# G(x) = x|x|:  G' = 2|x|,  G'' = 2 sign(x)
skew_deriv <- function(P) list(G = 2 * abs(P), G2 = 2 * sign(P))
skew_obj   <- function(L) colMeans(L * abs(L))

# Rank-r symmetric fastICA. Returns the n x r loadings t(U) %*% W plus the
# trace of sum(obj) over iterations, for convergence checking.
run_fastica <- function(U, r = nrow(U), n_iter = 50, seed = 1) {
  set.seed(seed)
  W     <- polar(matrix(rnorm(nrow(U) * r), nrow(U), r))
  trace <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    d <- skew_deriv(t(U) %*% W)
    W <- polar(U %*% d$G - sweep(W, 2, colSums(d$G2), "*"))
    trace[i] <- sum(skew_obj(t(U) %*% W))
  }
  list(L = t(U) %*% W, trace = trace)
}

# The objective is not concave and different random starts can land in
# different local optima, so we run several and keep the best. At K = 30 this
# matters: one start in five loses a real factor and ends up with a lower
# objective than the rest.
fit_one <- function(U_c, k) {
  nk     <- nrow(U_c)
  U_w    <- sqrt(nk) * t(U_c[, 1:k, drop = FALSE])   # whitened row space, k x nk
  starts <- lapply(seq_len(n_starts), function(s)
    run_fastica(U_w, r = k, n_iter = n_iter, seed = s))
  tot    <- sapply(starts, function(o) o$trace[n_iter])
  starts <- lapply(starts, function(o) {
    rownames(o$L) <- colnames(U_w)
    o
  })
  best <- starts[[which.max(tot)]]
  c(best, list(seed = which.max(tot), start_obj = tot,
               all_L = lapply(starts, `[[`, "L")))
}

fits <- lapply(K, function(k) fit_one(U_c, k))
names(fits) <- paste0("k", K)

# ---- individuals belonging to very sparse factors ----------------------------

# Rank-1 (unconstrained, unit-norm) update, all starts at once.
r1_update <- function(U, W) {
  P <- t(U) %*% W
  W <- U %*% (2 * abs(P)) - sweep(W, 2, colSums(2 * sign(P)), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

# Pool the K = 30 factors with the unconstrained maxima reached by running
# rank-1 from each of them. (The maxima found from 1000 random starts are a
# subset of these, so nothing is missed by leaving the random starts out here.)
U_w30 <- sqrt(n) * t(U_c[, 1:30, drop = FALSE])
L30   <- fits$k30$L
W_r1  <- t(U_c[, 1:30, drop = FALSE]) %*% L30 / sqrt(n)
for (i in 1:500) W_r1 <- r1_update(U_w30, W_r1)
pool <- cbind(L30, t(U_w30) %*% W_r1)

# A factor is "very sparse" if fewer than `sparse_nmax` individuals load above
# 2; an individual belongs to it if its loading exceeds `sparse_lmin`. Both
# choices are insensitive here: nmax in 20-25 and lmin in 5-8 give the same 40
# individuals. This picks out the relative pairs and the ITU/STU/GIH subgroups,
# but not the broader PJL subgroup (37 individuals above 2).
sparse_nmax <- 20
sparse_lmin <- 5

sparse_cols <- which(colSums(pool > 2) < sparse_nmax)
sparse_ind  <- which(apply(pool[, sparse_cols, drop = FALSE], 1, max) > sparse_lmin)
sparse_ids  <- rownames(Y)[sparse_ind]

# ---- rerun with those individuals removed ------------------------------------

keep_ind  <- setdiff(seq_len(n), sparse_ind)
no_sparse <- eigen_of(keep_ind)
rm(G_raw); gc()

fit_no_sparse <- fit_one(no_sparse$U_c, 30)

# ---- stability selection: split the SNPs by chromosome parity ----------------

# Column names are "chromosome:position", so we can split the SNPs into odd and
# even chromosomes and fit the two halves independently. Both use the pruned
# set of individuals, so the loadings are directly comparable.
chrom <- as.integer(sub(":.*", "", colnames(Y)))
stopifnot(!any(is.na(chrom)))

# Centred Gram matrix computed from scratch on a row/column submatrix.
eigen_of_block <- function(keep_rows, cols) {
  Yb <- Y[keep_rows, cols, drop = FALSE]
  nk <- nrow(Yb)
  mu <- colMeans(Yb)
  a  <- as.vector(Yb %*% mu)
  G  <- tcrossprod(Yb) -
          outer(a, rep(1, nk)) - outer(rep(1, nk), a) + sum(mu^2)
  G  <- (G + t(G)) / 2
  rm(Yb); gc()
  e <- eigen(G, symmetric = TRUE)
  U <- e$vectors
  rownames(U) <- rownames(Y)[keep_rows]
  list(d = sqrt(pmax(e$values, 0)), U_c = U[, 1:30, drop = FALSE])
}

halves <- list(odd  = which(chrom %% 2 == 1),
               even = which(chrom %% 2 == 0))
split_fits <- lapply(halves, function(cols) {
  e <- eigen_of_block(keep_ind, cols)
  list(n_snps = length(cols), d = e$d, U_c = e$U_c,
       fit = fit_one(e$U_c, 30))
})

rm(Y); gc()

saveRDS(list(fits      = fits,
             K         = K,
             d         = d,
             U_c       = U_c[, 1:max(K), drop = FALSE],
             n_iter    = n_iter,
             n_starts  = n_starts,
             no_sparse = list(dropped     = sparse_ids,
                              sparse_nmax = sparse_nmax,
                              sparse_lmin = sparse_lmin,
                              n_sparse_cols = length(sparse_cols),
                              d           = no_sparse$d,
                              U_c         = no_sparse$U_c[, 1:30, drop = FALSE],
                              fit         = fit_no_sparse),
             split     = split_fits),
        "output/tgp_fastica.rds")
