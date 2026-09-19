# Write a genotype matrix out in PLINK 1 binary format, so ADMIXTURE can read
# it. Originally written for data/tgp_data_matrix.rds; the input is now a
# parameter so the same code serves the HGDP data.
#
# Genotypes are counts of one allele (0, 1, 2), optionally with NA for missing.
# PLINK codes two bits per genotype, first sample in the low-order bits:
#   00 = hom A1, 01 = missing, 10 = het, 11 = hom A2
# so 0 -> 0, 1 -> 2, 2 -> 3, NA -> 1. Which allele is called A1 is irrelevant
# to ADMIXTURE's Q (swapping alleles just flips the corresponding row of P).
#
# Each SNP occupies ceil(n/4) bytes; when n is not a multiple of 4 the unused
# high-order bits of the last byte are zero-filled and ignored by the reader,
# which matters here because HGDP has 929 individuals.
#
# The input .rds may be either a plain matrix of genotypes, or the list written
# by code/prep_hgdp_data.R, whose $Y is a raw matrix of haplotypes (one row per
# chromosome, entries 0x30/0x31). With `pairs` set, consecutive rows are summed
# to give one 0/1/2 genotype per individual, which is what HGDP needs.
#
# An optional `thin=N` argument keeps every Nth SNP, which is how the
# ADMIXTURE runs are made tractable: ADMIXTURE's cost is linear in SNPs and
# tens of thousands of independent markers are ample to estimate ancestry
# proportions, whereas the fastICA fit wants all of them.
#
# Usage: Rscript code/write_plink_bed.R <input.rds> <outdir> <prefix> [pairs] [thin=N]

args   <- commandArgs(trailingOnly = TRUE)
inpath <- args[1]
outdir <- args[2]
prefix <- args[3]
pairs  <- "pairs" %in% args
thinarg <- grep("^thin=", args, value = TRUE)
thin   <- if (length(thinarg)) as.integer(sub("thin=", "", thinarg)) else 1L

x <- readRDS(inpath)
Y <- if (is.list(x)) x$Y else x
snps <- if (is.list(x) && !is.null(x$snps)) x$snps else NULL

if (thin > 1L) {
  keep_snps <- seq(1L, ncol(Y), by = thin)
  Y <- Y[, keep_snps, drop = FALSE]
  if (!is.null(snps)) snps <- snps[keep_snps, , drop = FALSE]
  cat("thinning to every", thin, "th SNP\n")
}

israw <- is.raw(Y[1, 1])
nrow_in <- nrow(Y); p <- ncol(Y)
n <- if (pairs) nrow_in / 2 else nrow_in
stopifnot(!pairs || nrow_in %% 2 == 0)

ids <- if (pairs) unique(sub("_[12]$", "", rownames(Y))) else rownames(Y)
stopifnot(length(ids) == n)

npad <- 4 * ceiling(n / 4)
cat(sprintf("writing %d individuals x %d SNPs (%d bytes/SNP)\n",
            n, p, npad / 4))

# ---- .fam and .bim ----------------------------------------------------------

write.table(data.frame(ids, ids, 0, 0, 0, -9),
            file.path(outdir, paste0(prefix, ".fam")), quote = FALSE,
            row.names = FALSE, col.names = FALSE)

snp_id <- colnames(Y)
bim <- data.frame(chr = sub("^chr", "", sub(":.*", "", snp_id)),
                  snp = snp_id, cm = 0,
                  bp  = as.numeric(sub(".*:", "", snp_id)),
                  a1  = if (is.null(snps)) "A" else snps$ref,
                  a2  = if (is.null(snps)) "G" else snps$alt)
write.table(bim, file.path(outdir, paste0(prefix, ".bim")), quote = FALSE,
            sep = "\t", row.names = FALSE, col.names = FALSE)

# ---- .bed -------------------------------------------------------------------

con <- file(file.path(outdir, paste0(prefix, ".bed")), "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)   # magic + SNP-major mode

chunk  <- 2000
n_miss <- 0
for (start in seq(1, p, by = chunk)) {
  cols <- start:min(start + chunk - 1, p)
  g    <- Y[, cols, drop = FALSE]
  if (israw) g <- matrix(as.integer(g) - 48L, nrow = nrow_in)
  if (pairs) g <- g[c(TRUE, FALSE), , drop = FALSE] +
                  g[c(FALSE, TRUE), , drop = FALSE]

  cc <- ifelse(is.na(g), 1L, ifelse(g == 0L, 0L, ifelse(g == 1L, 2L, 3L)))
  n_miss <- n_miss + sum(is.na(g))

  code <- matrix(0L, npad, length(cols))     # pad rows code as 00, ignored
  code[seq_len(n), ] <- cc
  dim(code) <- c(4L, npad / 4L, length(cols))
  bytes <- code[1, , , drop = FALSE] +  4L * code[2, , , drop = FALSE] +
      16L * code[3, , , drop = FALSE] + 64L * code[4, , , drop = FALSE]
  writeBin(as.raw(as.vector(bytes)), con)
  if (start %% 20000 == 1) cat(" ", start, "/", p, "\n")
}
close(con)
cat("missing genotypes:", n_miss, "\n")
cat("done:", file.info(file.path(outdir, paste0(prefix, ".bed")))$size / 1e6,
    "MB\n")
