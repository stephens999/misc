# log1p link Poisson NMF on the BBC news corpus, at c = 0.001.
#
# This is the model and the setting of Weine et al. (arXiv 2601.05845) --
# the c = 0.001 fit is the one their Figure shows -- run here on the same
# document-term matrix as code/fit_bbc_ica.R, so that the log1p NMF factors
# and the x|x| / EBproj factors can eventually be put side by side.
#
# Their fig_scripts/bbc_fig.R uses K = 10. That is kept, and K = 40 and 80 are
# added so there is a fit at each of the ranks the ICA analysis uses. The size
# factor s is the one code/fit_bbc_ica.R uses for log(1 + X/s): document
# length over mean document length.
#
# Nothing here feeds the write-up yet.
#
# Usage: BBC_NMF_K=10,40,80 Rscript code/fit_bbc_log1pnmf.R
# Reads  data/bbc_dtm.rds (from code/prep_bbc_data.R)
# Writes output/bbc_log1pnmf_c0.001_K<K>.rds, one per K.

suppressPackageStartupMessages({ library(Matrix); library(log1pNMF) })
if (!isTRUE(l10n_info()$`UTF-8`))
  for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8"))
    if (suppressWarnings(Sys.setlocale("LC_CTYPE", loc)) != "") break

cc  <- as.numeric(Sys.getenv("BBC_NMF_C", "0.001"))
Ks  <- as.integer(strsplit(Sys.getenv("BBC_NMF_K", "10,40,80"), ",")[[1]])

dat <- readRDS("data/bbc_dtm.rds")
Y   <- dat$X
s   <- Matrix::rowSums(Y); s <- s / mean(s)
cat(sprintf("%d docs x %d terms; c = %g; K = %s\n",
            nrow(Y), ncol(Y), cc, paste(Ks, collapse = ", ")))

for (K in Ks) {
  cat(sprintf("\n=== K = %d ===\n", K))
  set.seed(1)
  t0  <- Sys.time()
  fit <- fit_poisson_log1p_nmf(Y = Y, K = K, s = s, cc = cc,
                               loglik = "exact", init_method = "rank1")
  cat(sprintf("K = %d fitted in %.1f min\n", K,
              as.numeric(Sys.time() - t0, units = "mins")))
  colnames(fit$LL) <- colnames(fit$FF) <- paste0("k", seq_len(K))
  rownames(fit$LL) <- rownames(Y); rownames(fit$FF) <- colnames(Y)
  top <- apply(fit$FF, 2, function(f)
    paste(names(sort(f, decreasing = TRUE))[1:10], collapse = ", "))
  print(data.frame(factor = colnames(fit$FF), top_words = top), right = FALSE)
  saveRDS(list(cc = cc, K = K, fit = fit, s = s, labels = dat$labels,
               files = dat$files, terms = colnames(Y)),
          sprintf("output/bbc_log1pnmf_c%g_K%d.rds", cc, K))
  cat(sprintf("wrote output/bbc_log1pnmf_c%g_K%d.rds\n", cc, K))
}
