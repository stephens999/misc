# fastICA with the x|x| contrast, and EBproj with a fixed exponential prior, on
# the BBC news corpus of analysis/fastica_bbc.Rmd.
#
# The transformation is log(c + X/s): s is the document size factor (document
# length over mean document length, so order 1) when BBC_SCALE=1 and 1
# otherwise, and c is the pseudo-count, BBC_C, default 1. The columns are then
# centred and the data whitened to K principal components. Both contrasts run
# on the same whitening, so they are directly comparable.
#
# Note that log(c + z) = log(c) + log(1 + z/c), so the pseudo-count only
# rescales the data and shifts it by a constant. The constant is removed by
# the centring, so log(c + X/s) and log(1 + X/(cs)) give exactly the same
# whitening; the sparse form log1p(z/c) is what is computed. Small c pushes
# the transform towards presence/absence: at c = 0.001 a count of 1 maps to
# about 6.9 and a count of 2 to 7.6, so the gap between zero and one dwarfs
# the gaps between the non-zero counts.
#
# The search for local maxima extends the recipe of code/fit_immgent_ica.R in
# three ways:
#
#   1. The rank-r fit is run at r = K, five times from different random
#      starts, and *all* 5K columns are collected rather than only those of
#      the best fit. Duplicates among them are removed, and each distinct
#      rank-r direction is then relaxed with the rank-1 update -- the polar
#      step holds the columns of a rank-r fit mutually orthogonal, so as they
#      stand they are not maxima of the rank-1 objective. These are pooled
#      with 1000 random rank-1 starts, and all of the resulting rank-1 fits
#      are clustered together in one step.
#
#   2. The whole pipeline is run three times: on the full data and on each
#      half of a random split of the documents, so that a maximum can be
#      checked for reproducibility. The halves have their own whitenings, so
#      directions are compared in word space: a whitened-space direction w
#      corresponds to the word-space vector V[, 1:K] diag(d[1:K]) w, and the
#      cosine between two such vectors is comparable across fits.
#
#   3. Every solution set is then used to initialize every other one, in two
#      rounds. First across contrasts within a space, then across spaces
#      within a contrast -- a direction found in half 1 is carried into the
#      full data and into half 2 and relaxed there, and so on for all six
#      ordered pairs. A direction can be a genuine maximum in two places and
#      still be missed by one of the searches; this separates that case from
#      a direction that genuinely is not a maximum there. Anything new is
#      added to the set. Cross-space transfer goes through word space: from
#      the unit word-space vector g of the source, the nearest direction in
#      the target's row space is w = diag(1/d) V' g, renormalized.
#
# The data is small (2127 x 5861), so nothing here is batched or cached.
#
# Usage: BBC_K=40 BBC_SCALE=1 BBC_C=0.001 Rscript code/fit_bbc_ica.R
# Reads  data/bbc_dtm.rds (from code/prep_bbc_data.R)
# Writes output/bbc_ica_K<K>[_scaled][_c<c>].rds

suppressPackageStartupMessages({ library(Matrix) })
source("code/ebproj_pointexp.R")

# Term names carry a pound sign; keep the locale UTF-8 so it survives into the
# saved output. See the note in code/prep_bbc_data.R.
if (!isTRUE(l10n_info()$`UTF-8`))
  for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8"))
    if (suppressWarnings(Sys.setlocale("LC_CTYPE", loc)) != "") break

K        <- as.integer(Sys.getenv("BBC_K", "40"))
n_rankr  <- 5      # rank-r fits, each contributing K directions
n_iter   <- 500    # rank-r iterations
n_r1     <- 1000   # random rank-1 starts
n_r1iter <- 500    # rank-1 iterations
h_cut    <- 0.05   # clustering height: merge maxima correlated above 0.95
scaled   <- Sys.getenv("BBC_SCALE", "0") == "1"
pcount   <- as.numeric(Sys.getenv("BBC_C", "1"))   # pseudo-count c
# "docs": search for sparse directions over documents, the default and what
# the rest of the analysis uses. "words": search the row space instead, so a
# factor is a sparse set of words and the loadings are per word. The centring
# and whitening are identical; only which side of the SVD is searched
# changes. In word mode every fit lives in the same 5861-dimensional word
# space, so directions are directly comparable across the halves with no
# transfer needed.
space_of <- match.arg(Sys.getenv("BBC_SPACE", "docs"), c("docs", "words"))
# Which contrasts to run. The default pair is what the main analysis uses;
# adding "x^3" brings in the cubic skewness contrast. Every contrast is
# cross-initialized from every other, so the set matters.
cts_want <- strsplit(Sys.getenv("BBC_CONTRASTS", "x|x|,exp_fixed"), ",")[[1]]

out_file <- sprintf("output/bbc_ica_K%d%s%s%s%s.rds", K,
                    if (scaled) "_scaled" else "",
                    if (pcount == 1) "" else sprintf("_c%g", pcount),
                    if (space_of == "docs") "" else "_words",
                    if (identical(cts_want, c("x|x|", "exp_fixed"))) ""
                    else sprintf("_%dc", length(cts_want)))

dat <- readRDS("data/bbc_dtm.rds")
Y   <- dat$X

# Size factor, as in the log1p NMF analysis: document length relative to the
# average. Computed once on the full data, so that the halves stay on the same
# scale as each other and as the full fit.
sz <- Matrix::rowSums(Y); sz <- sz / mean(sz)
cat(sprintf("transform: log(%g + X%s); searching over %s\n",
            pcount, if (scaled) "/s" else "", space_of))

# ---- log1p, centre, whiten -------------------------------------------------
# log1p of a sparse matrix is sparse, so only the stored values change; the
# centring is kept out of the matrix and applied when the SVD is formed.

whiten <- function(rows) {
  L1 <- Y[rows, , drop = FALSE]
  if (scaled) L1 <- Diagonal(x = 1 / sz[rows]) %*% L1
  # log(c + z) up to the additive constant log(c), which the centring removes.
  L1@x <- log1p(L1@x / pcount)
  n  <- nrow(L1); p <- ncol(L1)
  cm <- Matrix::colMeans(L1); rm_ <- Matrix::rowMeans(L1)

  # Centre the side that is NOT being searched, so that the loadings the
  # contrasts see are mean zero. Searching over documents, the columns
  # (words) are centred: every column of Xc then sums to zero over documents,
  # so the column space is orthogonal to the all-ones vector in R^n and every
  # document-loading vector U'w is mean zero. Searching over words it is the
  # rows (documents) that must be centred, for the same reason in R^p.
  # Centring the wrong side leaves the all-ones direction in the space being
  # searched, and the contrasts will happily report it as a factor.
  Xc <- if (space_of == "docs") as.matrix(L1) - rep(cm, each = n)
        else                    as.matrix(L1) - rm_
  sv <- svd(Xc, nu = K, nv = K)
  list(n = n, rows = rows, d = sv$d, V = sv$v[, 1:K, drop = FALSE],
       pve = sv$d^2 / sum(Xc^2), col_mean = cm, row_mean = rm_,
       U_c = sv$u[, 1:K, drop = FALSE],
       L1 = L1, p = p,
       # K x m, rows scaled so that each has sum of squares m, where the
       # observations m are the documents in "docs" mode and the words in
       # "words" mode.
       U = if (space_of == "docs") {
             U <- sqrt(n) * t(sv$u[, 1:K, drop = FALSE])
             colnames(U) <- rownames(L1); U
           } else {
             U <- sqrt(p) * t(sv$v[, 1:K, drop = FALSE])
             colnames(U) <- colnames(L1); U
           })
}

# ---- the two contrasts, as interchangeable engines -------------------------
# Each supplies a rank-r fit, a rank-1 relaxation of an arbitrary set of
# starting directions, and the loadings of a direction. Everything downstream
# is shared.

polar <- function(W) { s <- svd(W); s$u %*% t(s$v) }

# G(x) = x|x|:  G' = 2|x|,  G'' = 2 sign(x)
xx_eng <- list(
  name = "x|x|",
  rankr = function(U, r, seed) {
    set.seed(seed)
    W <- polar(matrix(rnorm(nrow(U) * r), nrow(U), r))
    for (i in seq_len(n_iter)) {
      P <- crossprod(U, W)
      W <- polar(U %*% (2 * abs(P)) - sweep(W, 2, colSums(2 * sign(P)), "*"))
    }
    W
  },
  # Rank-1 (unconstrained, unit-norm) update applied to every start at once:
  # the columns never interact, so m runs cost one matrix product apiece.
  r1 = function(U, W) {
    for (i in seq_len(n_r1iter)) {
      W_old <- W
      P <- crossprod(U, W)
      W <- U %*% (2 * abs(P)) - sweep(W, 2, colSums(2 * sign(P)), "*")
      W <- sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
    }
    list(W = W, converged = 1 - abs(colSums(W * W_old)))
  },
  L   = function(U, W) crossprod(U, W),
  obj = function(U, W) colMeans({ L <- crossprod(U, W); L * abs(L) }))

# G(x) = x^3/3: G' = x^2, G'' = 2x. The classical skewness contrast, and the
# un-robustified sibling of x|x| = x^2 sign(x). Note that the fastICA Hessian
# correction, colSums(G''(P)) = 2 sum_i x_i, is essentially zero here: the
# data is centred on the side being searched, so the loadings are mean zero
# and the term vanishes. The update is therefore the plain fixed point
# W <- U P^2, renormalized.
x3_eng <- list(
  name = "x^3",
  rankr = function(U, r, seed) {
    set.seed(seed)
    W <- polar(matrix(rnorm(nrow(U) * r), nrow(U), r))
    for (i in seq_len(n_iter)) {
      P <- crossprod(U, W)
      W <- polar(U %*% (P^2) - sweep(W, 2, colSums(2 * P), "*"))
    }
    W
  },
  r1 = function(U, W) {
    for (i in seq_len(n_r1iter)) {
      W_old <- W
      P <- crossprod(U, W)
      W <- U %*% (P^2) - sweep(W, 2, colSums(2 * P), "*")
      W <- sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
    }
    list(W = W, converged = 1 - abs(colSums(W * W_old)))
  },
  L   = function(U, W) crossprod(U, W),
  obj = function(U, W) colMeans(crossprod(U, W)^3) / 3)

# EBproj, fixed exponential: pi0 = 0 and c^2 = sigma^2 = 1/2 held fixed. This
# is the `exp_fixed` configuration of code/fit_tgp_ebproj.R, the cheapest of
# the three and the one that worked best on HGDP and ImmGen-T. Loadings are
# the posterior means M(U'w), not U'w itself.
ef_par <- list(pi0 = 0, cc = sqrt(0.5), s2 = 0.5)
ef_eng <- list(
  name = "exp_fixed",
  rankr = function(U, r, seed) {
    set.seed(seed)
    pe_fit_r(U, r = r, n_iter = n_iter, pi0 = ef_par$pi0, cc = ef_par$cc,
             s2 = ef_par$s2, b = K/2, update_pi0 = FALSE, update_c = FALSE,
             update_s2 = FALSE)$W
  },
  r1 = function(U, W) {
    f <- pe_fit(U, W = W, n_iter = n_r1iter, tol = 1e-8, pi0 = ef_par$pi0,
                cc = ef_par$cc, s2 = ef_par$s2, b = K/2, update_pi0 = FALSE,
                update_c = FALSE, update_s2 = FALSE, g_method = "em")
    list(W = f$W, converged = f$last_chg, iters = f$iters,
         n_unconverged = f$n_unconverged)
  },
  L   = function(U, W) pe_stats(crossprod(U, W), ef_par$pi0, ef_par$cc, ef_par$s2)$M,
  obj = function(U, W) pe_objective(U, W, ef_par$pi0, ef_par$cc, ef_par$s2))

engines <- local({
  all <- list(xx_eng, ef_eng, x3_eng)
  names(all) <- sapply(all, `[[`, "name")
  stopifnot(cts_want %in% names(all))
  all[cts_want]
})
cat("contrasts:", paste(names(engines), collapse = ", "), "\n")

# ---- clustering ------------------------------------------------------------
# Maxima are merged when their loading vectors correlate above 1 - h. The
# correlation is *signed*: neither contrast is sign-symmetric, so w and -w are
# different solutions and must not be merged. For x|x| the loadings are U'w
# and U' has orthogonal columns scaled by sqrt(n), so this is exactly the
# inner product of the unit w's; for exp_fixed the posterior mean is a
# non-linear map of U'w and the two are not the same, so loadings it is.

cluster_L <- function(L, h = h_cut) {
  ok <- apply(L, 2, sd) > 1e-10
  hc <- hclust(as.dist(1 - cor(L[, ok, drop = FALSE])), method = "complete")
  cl_ok <- cutree(hc, h = h)
  cl <- rep(NA_integer_, ncol(L)); cl[ok] <- cl_ok
  sz <- sort(table(cl_ok), decreasing = TRUE)
  list(cl = cl, keep = which(ok)[match(as.integer(names(sz)), cl_ok)],
       sizes = as.integer(sz), ok = ok,
       sens = sapply(c(.90, .95, .99),
                     function(t) length(unique(cutree(hc, h = 1 - t)))))
}

# ---- the search ------------------------------------------------------------

find_maxima <- function(w, eng, verbose = TRUE) {
  U <- w$U
  say <- function(...) if (verbose) cat(sprintf(...))

  # (1) five rank-r fits, all K columns of each kept
  Wr <- do.call(cbind, lapply(seq_len(n_rankr), function(s) eng$rankr(U, K, s)))
  # The K columns within one fit are mutually orthogonal and so always
  # distinct; this removes directions repeated *between* fits.
  dr <- cluster_L(eng$L(U, Wr))
  Wr <- Wr[, dr$keep, drop = FALSE]
  say("  %d rank-r directions from %d fits, %d distinct\n",
      n_rankr * K, n_rankr, ncol(Wr))

  # (2) use each distinct rank-r direction to initialize the rank-1 update,
  # alongside 1000 random rank-1 starts. Both sets are run in the same call
  # and then clustered *jointly* below -- the rank-r-seeded runs are not a
  # separate pass, they are extra starts for the same search, and a maximum
  # reached from both routes must appear once, not twice.
  set.seed(1)
  W0 <- matrix(rnorm(K * n_r1), K, n_r1)
  W0 <- sweep(W0, 2, sqrt(colSums(W0^2)), "/")
  src <- c(rep("random", n_r1), rep("rank-r", ncol(Wr)))
  fit <- eng$r1(U, cbind(W0, Wr))
  say("  %d rank-1 runs, max 1 - |cos| at the last step %.2g\n",
      ncol(fit$W), max(fit$converged))

  # Joint clustering over all n_r1 + ncol(Wr) converged rank-1 fits.
  cl <- cluster_L(eng$L(U, fit$W))
  stopifnot(length(cl$cl) == ncol(fit$W), length(src) == ncol(fit$W))
  say("  %d maxima (at cuts .90/.95/.99: %s)\n",
      length(cl$keep), paste(cl$sens, collapse = "/"))

  list(W = fit$W[, cl$keep, drop = FALSE], src = src, cl = cl,
       sizes = cl$sizes, sens = cl$sens, n_rankr_distinct = ncol(Wr),
       n_random = sapply(cl$cl[cl$keep],
                         function(z) sum(cl$cl == z & src == "random", na.rm = TRUE)),
       n_from_rankr = sapply(cl$cl[cl$keep],
                             function(z) sum(cl$cl == z & src == "rank-r", na.rm = TRUE)),
       found_by = rep("search", length(cl$keep)),
       n_cross = list(), cross_cor = list())
}

# (3) cross-initialization: relax a set of directions found elsewhere -- by
# the other contrast, or in another space -- and see whether anything new
# turns up. `what` only names the source in the log.
cross_init <- function(w, eng, res, W_other, what = "other", verbose = TRUE) {
  U   <- w$U
  fit <- eng$r1(U, W_other)
  Lo  <- eng$L(U, res$W)
  Ln  <- eng$L(U, fit$W)
  ok  <- apply(Ln, 2, sd) > 1e-10
  bc  <- rep(0, ncol(Ln))
  bc[ok] <- apply(cor(Ln[, ok, drop = FALSE], Lo), 1, max)
  novel <- bc < 1 - h_cut
  if (verbose)
    cat(sprintf("  from %-18s %2d in: %2d land on an existing maximum, %d new\n",
                what, ncol(Ln), sum(!novel), sum(novel)))
  res$cross_cor[[what]] <- bc
  res$n_cross[[what]]   <- 0L
  if (!any(novel)) return(res)

  # Re-cluster the union so that several transferred directions landing on the
  # same new maximum are counted once. The new columns carry zero counts: they
  # were reached from the other contrast, not from a start of this one.
  n_own <- ncol(res$W)
  Wall  <- cbind(res$W, fit$W[, novel, drop = FALSE])
  cl    <- cluster_L(eng$L(U, Wall))
  add   <- cl$keep[cl$keep > n_own]
  res$W            <- cbind(res$W, Wall[, add, drop = FALSE])
  res$sizes        <- c(res$sizes, rep(0L, length(add)))
  res$n_random     <- c(res$n_random, rep(0L, length(add)))
  res$n_from_rankr <- c(res$n_from_rankr, rep(0L, length(add)))
  res$found_by     <- c(res$found_by, rep(what, length(add)))
  res$n_cross[[what]] <- length(add)
  res
}

# Attach the derived quantities: loadings, objective, the word-space
# representative used for comparison across spaces, and the word loadings.
finalize <- function(w, eng, o) {
  o$L   <- eng$L(w$U, o$W); rownames(o$L) <- colnames(w$U)
  o$obj <- eng$obj(w$U, o$W)
  # Word-space representative, used to compare directions across fits with
  # different whitenings. In word mode the loadings are already word-space
  # vectors, so no reconstruction is needed.
  G <- if (space_of == "docs") w$V %*% (w$d[1:K] * o$W) else o$L
  o$G <- sweep(G, 2, sqrt(colSums(G^2)), "/")
  rownames(o$G) <- colnames(Y)
  o$words <- cross_loadings(w, o$L)
  o
}

# Carry directions from one space into another through word space: the unit
# word-space vector g of the source maps to the nearest direction in the
# target's row space, w = diag(1/d) V' g, renormalized.
transfer <- function(w, G) {
  # docs mode: g = V diag(d) w, so w = diag(1/d) V' g.
  # words mode: the loadings are sqrt(p) V w, so w = V' g up to scale.
  W <- if (space_of == "docs") (1 / w$d[1:K]) * crossprod(w$V, G)
       else crossprod(w$V, G)
  sweep(W, 2, sqrt(colSums(W^2)) + 1e-15, "/")
}

# Stage 1: the search, independently in each space and under each contrast.
search_all <- function(spaces) {
  out <- list()
  for (nm in names(spaces)) {
    cat(sprintf("\n=== %s: %d documents, K = %d ===\n", nm, spaces[[nm]]$n, K))
    for (eng in engines) {
      cat(sprintf(" %s\n", eng$name))
      out[[nm]][[eng$name]] <- finalize(spaces[[nm]], eng,
                                        find_maxima(spaces[[nm]], eng))
    }
  }
  out
}

# Stage 2: cross-contrast initialization, within each space.
cross_contrasts <- function(spaces, out) {
  cat("\n=== cross-contrast initialization ===\n")
  for (nm in names(spaces)) {
    cat(sprintf(" %s\n", nm))
    W <- lapply(out[[nm]], `[[`, "W")
    for (eng in engines) {
      o <- out[[nm]][[eng$name]]
      for (other in setdiff(names(W), eng$name))
        o <- cross_init(spaces[[nm]], eng, o, W[[other]], other)
      out[[nm]][[eng$name]] <- finalize(spaces[[nm]], eng, o)
    }
  }
  out
}

# Stage 3: cross-space initialization, within each contrast. Every space is
# started from every other space's maxima. The sources are the sets as they
# stand after stage 2, so the round is symmetric and order does not matter.
cross_spaces <- function(spaces, out) {
  cat("\n=== cross-space initialization ===\n")
  G0 <- lapply(out, function(s) lapply(s, `[[`, "G"))
  for (nm in names(spaces)) {
    cat(sprintf(" %s\n", nm))
    for (eng in engines) {
      o <- out[[nm]][[eng$name]]
      for (src in setdiff(names(spaces), nm))
        o <- cross_init(spaces[[nm]], eng, o,
                        transfer(spaces[[nm]], G0[[src]][[eng$name]]), src)
      out[[nm]][[eng$name]] <- finalize(spaces[[nm]], eng, o)
    }
  }
  out
}

# ---- word loadings ---------------------------------------------------------
# Regress the centred log1p matrix on the non-negative part of each factor's
# document loadings, scaled to a maximum of 1, plus an intercept. This is what
# code/fit_immgent_ebproj.R does, and it puts the words that go with a factor
# on the same scale as the data.

# Regress the transformed matrix on the non-negative part of the loadings
# (scaled to a maximum of 1), with an intercept. In docs mode the loadings are
# per document and this returns word loadings; in words mode the loadings are
# per word and it returns document loadings. Same computation either way, on
# the matrix oriented so that the loadings index its rows.
cross_loadings <- function(w, L) {
  M  <- if (space_of == "docs") w$L1 else Matrix::t(w$L1)
  nr <- nrow(M)
  ctr <- if (space_of == "docs") w$col_mean else w$row_mean
  A <- pmax(L, 0); mx <- apply(A, 2, max); A <- sweep(A, 2, mx, "/")
  G <- rbind(cbind(nr, t(colSums(A))), cbind(colSums(A), crossprod(A)))
  B <- solve(G, rbind(nr * ctr, as.matrix(crossprod(A, M))))
  colnames(B) <- colnames(M)
  list(B = B[-1, , drop = FALSE], B0 = B[1, ], scale = mx,
       above_0.5 = colSums(A > 0.5))
}

# ---- run on the full data and on each half ---------------------------------

set.seed(1)
half <- sample(rep(1:2, length.out = nrow(Y)))
spaces <- list(full  = whiten(seq_len(nrow(Y))),
               half1 = whiten(which(half == 1)),
               half2 = whiten(which(half == 2)))
for (nm in names(spaces))
  cat(sprintf("%s: %d docs, top %d PCs explain %.1f%% of the variance\n",
              nm, spaces[[nm]]$n, K, 100 * sum(spaces[[nm]]$pve[1:K])))

t0  <- Sys.time()
res <- search_all(spaces)
res <- cross_contrasts(spaces, res)
res <- cross_spaces(spaces, res)
cat("\n=== maxima per space and contrast ===\n")
print(sapply(res, function(s) sapply(s, function(f) ncol(f$W))))
cat(sprintf("total %.1f min\n", as.numeric(Sys.time() - t0, units = "mins")))

saveRDS(list(K = K, n = nrow(Y), p = ncol(Y), h_cut = h_cut,
             n_rankr = n_rankr, n_r1 = n_r1, n_r1iter = n_r1iter,
             n_iter = n_iter, half = half,
             scaled = scaled, pcount = pcount, space = space_of,
             contrasts = names(engines),
             labels = dat$labels, files = dat$files, size = sz,
             terms = colnames(Y),
             spaces = lapply(spaces, function(w)
               w[c("n", "rows", "d", "V", "pve", "col_mean", "U_c")]),
             fits = res),
        out_file)
cat("wrote", out_file, "\n")
