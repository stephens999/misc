# KING-robust kinship among the 929 HGDP individuals, from the same SNPs the
# fastICA analysis uses. This is here to answer two questions the analysis
# raises: whether the HGDP WGS panel still contains close relatives, and
# whether the individuals that the x|x| contrast singles out into very narrow
# factors are related or are genuine fine-scale population substructure.
#
# The KING-robust estimator for a pair (i, j) is
#   phi_ij = (N_AaAa - 2 N_AAaa) / (N_Aa_i + N_Aa_j)
# with the usual thresholds: > 0.354 duplicate/MZ, > 0.177 first-degree,
# > 0.0884 second-degree, > 0.0442 third-degree. It is computed from the
# haplotype matrix, summing chromosome pairs to genotypes, in SNP blocks.
#
# Caveat worth remembering when reading the output: KING-robust assumes the
# pair shares ancestry, and in small strongly drifted populations (Surui,
# Karitiana, Pima, Bougainville) background relatedness inflates it, so a
# "second-degree" call there need not be a pedigree relationship.
#
# Usage: Rscript code/fit_hgdp_kinship.R
# Writes output/hgdp_kinship.rds.

meta <- readRDS("data/hgdp_meta.rds")
dat  <- readRDS("data/hgdp_hap_matrix.rds")
Y    <- dat$Y
nh <- nrow(Y); p <- ncol(Y); n <- nh / 2
ids <- unique(sub("_[12]$", "", rownames(Y)))
stopifnot(identical(ids, as.character(meta$sample)))

blocks <- split(seq_len(p), ceiling(seq_len(p) / 20000))
Nhh <- matrix(0, n, n); Nopp <- matrix(0, n, n); nhet <- numeric(n)
for (cols in blocks) {
  h <- matrix(as.integer(Y[, cols, drop = FALSE]) - 48L, nrow = nh)
  g <- h[c(TRUE, FALSE), , drop = FALSE] + h[c(FALSE, TRUE), , drop = FALSE]
  H <- (g == 1) * 1.0; A0 <- (g == 0) * 1.0; A2 <- (g == 2) * 1.0
  Nhh  <- Nhh + tcrossprod(H)
  Nopp <- Nopp + tcrossprod(A0, A2) + tcrossprod(A2, A0)
  nhet <- nhet + rowSums(H)
}
phi <- (Nhh - 2 * Nopp) / outer(nhet, nhet, "+")
dimnames(phi) <- list(ids, ids)
diag(phi) <- NA

deg <- function(x)
  ifelse(x > 0.354, "dup/MZ", ifelse(x > 0.177, "1st",
  ifelse(x > 0.0884, "2nd", ifelse(x > 0.0442, "3rd", "unrelated"))))

ut <- which(upper.tri(phi), arr.ind = TRUE)
pairs <- data.frame(a = ids[ut[, 1]], b = ids[ut[, 2]],
                    pop_a = as.character(meta$pop)[ut[, 1]],
                    pop_b = as.character(meta$pop)[ut[, 2]],
                    kinship = phi[ut])
pairs <- pairs[pairs$kinship > 0.0442, ]
pairs$degree <- deg(pairs$kinship)
pairs <- pairs[order(-pairs$kinship), ]
rownames(pairs) <- NULL

counts <- table(factor(deg(phi[ut]),
                       levels = c("dup/MZ", "1st", "2nd", "3rd", "unrelated")))
print(counts)
cat("individuals in a pair at 2nd degree or closer:",
    length(unique(c(pairs$a, pairs$b)[pairs$kinship > 0.0884])), "\n")

saveRDS(list(phi = phi, ids = ids, max_kinship = apply(phi, 1, max, na.rm = TRUE),
             pairs = pairs, counts = counts, n_pairs = nrow(ut)),
        "output/hgdp_kinship.rds")
cat("wrote output/hgdp_kinship.rds\n")
