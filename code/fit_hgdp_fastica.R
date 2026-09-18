# fastICA with the x|x| contrast on the HGDP phased haplotypes, following the
# 1000 Genomes analysis in analysis/fastica_1kg.Rmd and code/fit_tgp_fastica.R.
#
# The data are 2 x 929 = 1858 haplotypes by ~1.7e5 thinned autosomal SNPs,
# coded as derived-allele indicators (data/hgdp_hap_matrix.rds, built by
# code/prep_hgdp_data.R). Columns are centred and nothing else is done to them.
#
# What this does, in the order the analysis reads:
#   1. centred SVD of the haplotype matrix, via the 1858 x 1858 Gram matrix;
#   2. for each K in `K`: a rank-r fastICA from five random starts, best
#      objective kept, and 1000 independent rank-1 runs from random starts in
#      the K-dimensional whitened space, clustered into distinct local maxima;
#   3. a stability check: the SNPs are split by chromosome parity and the
#      K = 40 fit and rank-1 search redone independently on each half;
#   4. a check on whether pruning related/outlier individuals is needed. On
#      these data it is not -- see analysis/fastica_hgdp.Rmd -- so the fits
#      above use every individual, and the pruned rerun is kept only as the
#      evidence for that claim. It is done at K = 40 only.
#
# Usage: Rscript code/fit_hgdp_fastica.R
# Writes output/hgdp_fastica.rds.

K        <- c(40, 60)   # whitening dimensions
n_iter   <- 200    # rank-r iterations
n_starts <- 5      # random starts for the rank-r fit
n_r1     <- 1000   # rank-1 runs
n_r1iter <- 500    # rank-1 iterations

meta        <- readRDS("data/hgdp_meta.rds")
meta_pop    <- as.character(meta$pop)
meta_sample <- meta$sample

dat <- readRDS("data/hgdp_hap_matrix.rds")
Y   <- dat$Y                       # raw, 1858 x nsnp, entries 0x30 / 0x31
n   <- nrow(Y)
p   <- ncol(Y)
cat(sprintf("%d haplotypes x %d SNPs\n", n, p))

# Column blocks, so that no more than `blk` SNPs are numeric at a time.
blk    <- 20000
blocks <- split(seq_len(p), ceiling(seq_len(p) / blk))

as_num <- function(cols) {
  matrix(as.numeric(Y[, cols, drop = FALSE]) - 48, nrow = n)
}

# ---- centred SVD via the Gram matrix ----------------------------------------
#
# As in the TGP script we never form a centred copy of Y. With mu the column
# means over the rows we are keeping and a = Y[keep, ] mu,
#   C C' = (YY')[keep, keep] - a 1' - 1 a' + (mu'mu) 11'.
# Here YY' and the column sums are accumulated over SNP blocks.

G_raw <- matrix(0, n, n)
csums <- numeric(p)
for (cols in blocks) {
  Yb          <- as_num(cols)
  G_raw       <- G_raw + tcrossprod(Yb)
  csums[cols] <- colSums(Yb)
}
G_raw <- (G_raw + t(G_raw)) / 2

# a = Y[keep, ] mu, accumulated the same way.
a_of <- function(keep, mu) {
  a <- numeric(length(keep))
  for (cols in blocks) a <- a + as_num(cols)[keep, , drop = FALSE] %*% mu[cols]
  as.vector(a)
}

eigen_of <- function(keep) {
  nk <- length(keep)
  mu <- if (nk == n) csums / n else {
    drop <- setdiff(seq_len(n), keep)
    cd   <- numeric(p)
    for (cols in blocks) cd[cols] <- colSums(as_num(cols)[drop, , drop = FALSE])
    (csums - cd) / nk
  }
  a <- a_of(keep, mu)
  G <- G_raw[keep, keep, drop = FALSE] -
         outer(a, rep(1, nk)) - outer(rep(1, nk), a) + sum(mu^2)
  G <- (G + t(G)) / 2
  e <- eigen(G, symmetric = TRUE)
  U <- e$vectors
  rownames(U) <- rownames(Y)[keep]
  list(d = sqrt(pmax(e$values, 0)), U_c = U)
}

full <- eigen_of(seq_len(n))

# ---- fastICA ----------------------------------------------------------------

polar <- function(W) {
  eig <- eigen(t(W) %*% W, symmetric = TRUE)
  W %*% (eig$vectors %*% diag(1 / sqrt(pmax(eig$values, 1e-14))) %*% t(eig$vectors))
}

# G(x) = x|x|:  G' = 2|x|,  G'' = 2 sign(x)
skew_deriv <- function(P) list(G = 2 * abs(P), G2 = 2 * sign(P))
skew_obj   <- function(L) colMeans(L * abs(L))

run_fastica <- function(U, r = nrow(U), n_iter = 200, seed = 1) {
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

whiten <- function(U_c, k) sqrt(nrow(U_c)) * t(U_c[, 1:k, drop = FALSE])

fit_one <- function(U_c, k) {
  U_w    <- whiten(U_c, k)
  starts <- lapply(seq_len(n_starts), function(s)
    run_fastica(U_w, r = k, n_iter = n_iter, seed = s))
  starts <- lapply(starts, function(o) { rownames(o$L) <- colnames(U_w); o })
  tot    <- sapply(starts, function(o) o$trace[n_iter])
  c(starts[[which.max(tot)]],
    list(seed = which.max(tot), start_obj = tot,
         all_L = lapply(starts, `[[`, "L")))
}

# Rank-1 (unconstrained, unit-norm) update applied to every start at once: the
# columns of W never interact, so 1000 runs cost one matrix product apiece.
r1_update <- function(U, W) {
  P <- t(U) %*% W
  W <- U %*% (2 * abs(P)) - sweep(W, 2, colSums(2 * sign(P)), "*")
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

# 1000 random rank-1 runs, plus rank-1 started at each rank-r factor; the two
# sets are clustered together so a maximum found by either route appears once.
maxima_of <- function(U_c, k, L_rankr, seed = 1) {
  nk  <- nrow(U_c)
  U_w <- whiten(U_c, k)

  set.seed(seed)
  W <- matrix(rnorm(k * n_r1), k, n_r1)
  W <- sweep(W, 2, sqrt(colSums(W^2)), "/")
  for (i in seq_len(n_r1iter)) { W_old <- W; W <- r1_update(U_w, W) }
  converged <- max(1 - abs(colSums(W * W_old)))

  W_k <- t(U_c[, 1:k, drop = FALSE]) %*% L_rankr / sqrt(nk)
  for (i in seq_len(n_r1iter)) W_k <- r1_update(U_w, W_k)

  L  <- cbind(t(U_w) %*% W, t(U_w) %*% W_k)
  # Signed correlation: x|x| is not sign-symmetric, so w and -w are different
  # solutions and must not be merged.
  cl <- cutree(hclust(as.dist(1 - cor(L)), method = "complete"), h = 0.05)
  j  <- match(unique(cl), cl)
  M  <- L[, j, drop = FALSE]
  rownames(M) <- rownames(U_c)
  list(L          = M,
       n_random   = sapply(unique(cl), function(z) sum(cl[seq_len(n_r1)] == z)),
       from_rankr = sapply(unique(cl), function(z) any(cl[-seq_len(n_r1)] == z)),
       converged  = converged,
       n_cut      = sapply(c(0.90, 0.95, 0.99), function(tau)
                       length(unique(cutree(hclust(as.dist(1 - cor(L)),
                                                   method = "complete"),
                                            h = 1 - tau)))))
}

fits <- lapply(K, function(k) {
  cat(sprintf("K = %d: rank-r fit\n", k))
  fk <- fit_one(full$U_c, k)
  cat(sprintf("K = %d: maxima\n", k))
  mk <- maxima_of(full$U_c, k, fk$L)
  cat(sprintf("K = %d: %d maxima (%d reached only from the rank-r factors)\n",
              k, ncol(mk$L), sum(mk$n_random == 0)))
  list(K = k, fit = fk, maxima = mk)
})
names(fits) <- paste0("k", K)

# ---- stability: split the SNPs by chromosome parity -------------------------
#
# The two halves are independent sets of markers, so a maximum that both find
# is not an artifact of which SNPs were used. Column names are "chrN:pos", so
# splitting on chromosome parity keeps whole chromosomes together and cannot
# split residual LD across the two halves.

chrom <- as.integer(sub("^chr", "", sub(":.*", "", colnames(Y))))
stopifnot(!any(is.na(chrom)))
halves <- list(odd = which(chrom %% 2 == 1), even = which(chrom %% 2 == 0))

# Centred Gram matrix over an arbitrary set of SNP columns, accumulated in
# blocks exactly as for the full data.
gram_cols <- function(cols) {
  cb <- split(cols, ceiling(seq_along(cols) / blk))
  G  <- matrix(0, n, n)
  mu <- vector("list", length(cb))
  for (i in seq_along(cb)) {
    Yb     <- matrix(as.numeric(Y[, cb[[i]], drop = FALSE]) - 48, nrow = n)
    G      <- G + tcrossprod(Yb)
    mu[[i]] <- colMeans(Yb)
  }
  a <- numeric(n)
  for (i in seq_along(cb)) {
    Yb <- matrix(as.numeric(Y[, cb[[i]], drop = FALSE]) - 48, nrow = n)
    a  <- a + as.vector(Yb %*% mu[[i]])
  }
  s2 <- sum(unlist(mu)^2)
  G  <- G - outer(a, rep(1, n)) - outer(rep(1, n), a) + s2
  (G + t(G)) / 2
}

K_split <- 40
split_fits <- lapply(names(halves), function(h) {
  cols <- halves[[h]]
  cat(sprintf("%s chromosomes: %d SNPs\n", h, length(cols)))
  e <- eigen(gram_cols(cols), symmetric = TRUE)
  U <- e$vectors
  rownames(U) <- rownames(Y)
  fk <- fit_one(U, K_split)
  mk <- maxima_of(U, K_split, fk$L)
  cat(sprintf("%s: %d maxima (%d only from the rank-r factors)\n",
              h, ncol(mk$L), sum(mk$n_random == 0)))
  # all_L is five 1858 x 40 matrices per fit and is not used downstream.
  fk$all_L <- NULL
  list(n_snps = length(cols), d = sqrt(pmax(e$values, 0)),
       U_c = U[, 1:K_split, drop = FALSE], fit = fk, maxima = mk)
})
names(split_fits) <- names(halves)

# ---- is pruning needed? -----------------------------------------------------
#
# The 1000 Genomes analysis found that a handful of related individuals and
# sample outliers distorted the whole landscape, and pruned them. Here that
# turns out not to be necessary, and this block is the evidence: it runs the
# screen and the pruned rerun at K = 40 so the analysis can show that removing
# the individuals changes nothing but the two within-population maxima.
#
# The 1000 Genomes rule -- "fewer than 20 individuals load above 2" -- does not
# transfer, because every 1000 Genomes population has 60-100 individuals while
# HGDP populations run from 6 to 46. Used as written it deletes whole
# populations: San (6), Colombian (7), Surui (8) and Mbuti (13) each have their
# own maximum and each looks "very sparse" by that rule. So the test is made
# relative to population size: a maximum is screened in if it is narrow, and
# judged a nuisance only if the individuals it singles out are a minority of
# the population most of them belong to.

K_prune        <- 40
sparse_nmax    <- 40     # haplotypes loading above 2: the factor must be narrow
sparse_lmin    <- 5      # loading at which a haplotype "belongs to" a factor
sparse_fracmax <- 0.5    # ... and must be a minority of its own population

ind_of <- sub("_[12]$", "", rownames(Y))
pop_of <- meta_pop[match(ind_of, meta_sample)]

pool   <- cbind(fits[[paste0("k", K_prune)]]$fit$L,
                fits[[paste0("k", K_prune)]]$maxima$L)
narrow <- which(colSums(pool > 2) < sparse_nmax)

# For each narrow maximum: the individuals it picks out, the population most of
# them come from, and what fraction of that population they are.
screen <- lapply(narrow, function(j) {
  who <- unique(ind_of[pool[, j] > sparse_lmin])
  if (length(who) == 0) return(NULL)
  tp  <- sort(table(pop_of[match(who, ind_of)]), decreasing = TRUE)
  top <- names(tp)[1]
  # meta_pop is per individual; pop_of is per haplotype, so the denominator
  # has to come from meta_pop or every fraction is halved.
  n_pop <- sum(meta_pop == top)
  list(col = j, who = who, top_pop = top, n_in = tp[[1]],
       n_pop = n_pop, frac = tp[[1]] / n_pop)
})
screen <- Filter(Negate(is.null), screen)
stopifnot(length(screen) > 0)

screen_tab <- do.call(rbind, lapply(screen, function(z)
  data.frame(n_individuals = length(z$who), top_pop = z$top_pop,
             n_in_pop = z$n_in, pop_size = z$n_pop, frac = round(z$frac, 2),
             obj = round(skew_obj(pool[, z$col, drop = FALSE]), 2),
             n_above2 = sum(pool[, z$col] > 2))))

is_sparse  <- sapply(screen, function(z) z$frac < sparse_fracmax)
sparse_ids <- unique(unlist(lapply(screen[is_sparse], `[[`, "who")))
keep       <- which(!(ind_of %in% sparse_ids))
cat(sprintf("prune check: %d narrow maxima, %d within-population, %d individuals, %d haplotypes left\n",
            nrow(screen_tab), sum(is_sparse), length(sparse_ids), length(keep)))

pruned      <- eigen_of(keep)
fit_pruned  <- fit_one(pruned$U_c, K_prune)
max_pruned  <- maxima_of(pruned$U_c, K_prune, fit_pruned$L)
cat(sprintf("prune check: %d maxima after pruning\n", ncol(max_pruned$L)))

saveRDS(list(
  K        = K,
  n_iter   = n_iter,
  n_starts = n_starts,
  n_r1     = n_r1,
  n_r1iter = n_r1iter,
  nsnp     = p,
  hap_ids  = rownames(Y),
  d        = full$d,
  U_c      = full$U_c[, 1:max(K), drop = FALSE],
  fits     = fits,
  split    = c(split_fits, list(K = K_split)),
  prune_check = list(K              = K_prune,
                     sparse_nmax    = sparse_nmax,
                     sparse_lmin    = sparse_lmin,
                     sparse_fracmax = sparse_fracmax,
                     screen_tab     = screen_tab,
                     dropped        = sparse_ids,
                     hap_ids        = rownames(Y)[keep],
                     maxima         = max_pruned)),
  "output/hgdp_fastica.rds")
cat("wrote output/hgdp_fastica.rds\n")
