# Does the EBproj solution make a good starting point for a topic model?
#
# Takes the normalized loadings of the fixed-exponential EBproj fit at K = 80
# on the full data -- the memberships used in the structure plots, i.e. the
# positive part of the posterior means scaled to a maximum of 1 -- and uses
# them to initialize a fastTopics Poisson NMF on the raw counts. The same fit
# is then run from fastTopics' own default initialization at the same K, so
# the two can be compared on likelihood and on which topics they find.
#
# The point is that EBproj works on a centred, whitened, log1p-transformed
# matrix and produces signed factors, while the topic model works on the raw
# counts and produces non-negative ones. If the EBproj maxima are picking out
# real structure, they should be a good place to start the topic model from.
#
# Not part of the write-up yet.
#
# Usage: Rscript code/fit_bbc_fasttopics.R
# Reads  data/bbc_dtm.rds, output/bbc_ica_K80_scaled.rds
# Writes output/bbc_fasttopics_K80.rds

suppressPackageStartupMessages({ library(Matrix); library(fastTopics) })
if (!isTRUE(l10n_info()$`UTF-8`))
  for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8"))
    if (suppressWarnings(Sys.setlocale("LC_CTYPE", loc)) != "") break

numiter <- as.integer(Sys.getenv("BBC_FT_ITER", "500"))
src     <- Sys.getenv("BBC_FT_SRC", "output/bbc_ica_K80_scaled.rds")
cfg     <- Sys.getenv("BBC_FT_CFG", "exp_fixed")

dat <- readRDS("data/bbc_dtm.rds")
X   <- dat$X
r   <- readRDS(src)
ef  <- r$fits$full[[cfg]]

# Normalized EBproj memberships: positive part, each column scaled to max 1.
A <- pmax(ef$L, 0)
A <- sweep(A, 2, apply(A, 2, max), "/")
K <- ncol(A)
cat(sprintf("%d docs x %d terms; init from %s %s, K = %d\n",
            nrow(X), ncol(X), basename(src), cfg, K))
cat(sprintf("memberships: %.1f%% exactly zero; docs above 0.5 per factor, median %.0f\n",
            100 * mean(A == 0), median(colSums(A > 0.5))))

# fastTopics needs strictly positive loadings to start from, and a document
# that loads on nothing would be stuck at zero, so floor at a small value.
L0 <- pmax(A, 1e-4)

# With L supplied and F left to be filled in, fastTopics requires the random
# initializer for F; topicscore is only available when neither is given.
set.seed(1)
t0 <- Sys.time()
fit0_ebp <- init_poisson_nmf(X, L = L0, init.method = "random", verbose = "none")
ebp <- fit_poisson_nmf(X, fit0 = fit0_ebp, numiter = numiter, method = "scd",
                       verbose = "none")
cat(sprintf("EBproj-initialized: %.1f min, loglik %.1f\n",
            as.numeric(Sys.time() - t0, units = "mins"), sum(loglik_poisson_nmf(X, ebp))))

set.seed(1)
t0 <- Sys.time()
def <- fit_poisson_nmf(X, k = K, numiter = numiter, method = "scd",
                       init.method = "topicscore", verbose = "none")
cat(sprintf("default init:       %.1f min, loglik %.1f\n",
            as.numeric(Sys.time() - t0, units = "mins"), sum(loglik_poisson_nmf(X, def))))

# ---- compare ---------------------------------------------------------------

top_words <- function(F, j, n = 10)
  paste(rownames(F)[order(-F[, j])[1:n]], collapse = ", ")

for (f in list(ebp, def)) { rownames(f$F) <- colnames(X); rownames(f$L) <- rownames(X) }
rownames(ebp$F) <- rownames(def$F) <- colnames(X)
rownames(ebp$L) <- rownames(def$L) <- rownames(X)

# Mutual nearest neighbours between the two fits, on the topic-word
# distributions, which is the comparison that does not depend on scaling.
mnn <- function(C) {
  a <- max.col(C, "first"); b <- max.col(t(C), "first"); i <- seq_len(nrow(C))
  data.frame(i = i, j = a, cos = C[cbind(i, a)], mutual = b[a] == i)
}
pn <- function(M) sweep(M, 2, sqrt(colSums(M^2)), "/")
C  <- crossprod(pn(ebp$F), pn(def$F))
m  <- mnn(C)
cat(sprintf("\n%d of %d EBproj-initialized topics are mutual nearest neighbours of a default topic; %d above 0.7\n",
            sum(m$mutual), K, sum(m$mutual & m$cos > 0.7)))
print(data.frame(
  ebproj_init = sapply(seq_len(K), function(j) top_words(ebp$F, j, 6)),
  default     = ifelse(m$mutual, sapply(m$j, function(j) top_words(def$F, j, 6)), "-"),
  cos         = ifelse(m$mutual, round(m$cos, 2), NA)), right = FALSE)

# How much of the initialization survived? Not a column-by-column question:
# F starts random, so the topic labels shuffle in the first few updates even
# though the loading structure carries over. Match the starting memberships to
# the fitted loadings the same way as everything else, and do it against the
# default fit too as a baseline -- if the two agree equally well with the
# starting point, the initialization left no trace.
init_match <- function(fit, what) {
  m <- mnn(crossprod(pn(A), pn(fit$L)))
  cat(sprintf("init -> %-8s mutual %2d of %2d, cosine median %.2f (%.2f-%.2f)\n",
              what, sum(m$mutual), K, median(m$cos[m$mutual]),
              min(m$cos[m$mutual]), max(m$cos[m$mutual])))
  m
}
cat("\n")
mi_ebp <- init_match(ebp, "ebproj")
mi_def <- init_match(def, "default")

saveRDS(list(K = K, numiter = numiter, src = src, cfg = cfg,
             ebp = ebp, def = def, A_init = A,
             loglik_ebp = sum(loglik_poisson_nmf(X, ebp)),
             loglik_def = sum(loglik_poisson_nmf(X, def)),
             match = m, init_match_ebp = mi_ebp, init_match_def = mi_def,
             labels = dat$labels, terms = colnames(X)),
        "output/bbc_fasttopics_K80.rds")
cat("wrote output/bbc_fasttopics_K80.rds\n")
