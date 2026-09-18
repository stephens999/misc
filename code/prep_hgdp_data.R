# Build the HGDP metadata and haplotype matrix used by analysis/fastica_hgdp.Rmd.
#
# Inputs are produced outside R:
#   * the sample metadata table from
#     https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/metadata/
#   * the statistically phased autosomal VCF from
#     https://ngs.sanger.ac.uk/production/hgdp/hgdp_wgs.20190516/statphase/,
#     reduced by code/vcf_thin.c to a thinned set of biallelic SNPs written as
#     derived-allele haplotype calls (see that file for the filters):
#
#       cc -O2 -o vcf_thin code/vcf_thin.c
#       gzip -dc hgdp_wgs.20190516.statphase.autosomes.vcf.gz \
#         | ./vcf_thin 6000 0.05 hgdp_thin
#
#     6000 is the minimum bp between kept SNPs and gives 215468 SNPs; 12000
#     gives 140628, and analysis/fastica_hgdp.Rmd shows the two agree.
#
# Usage: Rscript code/prep_hgdp_data.R <thin_prefix> <meta_txt>
# Writes data/hgdp_meta.rds and data/hgdp_hap_matrix.rds.

args   <- commandArgs(trailingOnly = TRUE)
prefix <- if (length(args) > 0) args[1] else "hgdp_thin"
metaf  <- if (length(args) > 1) args[2] else "hgdp_wgs.20190516.metadata.txt"

# ---- metadata ---------------------------------------------------------------

m <- read.delim(metaf, stringsAsFactors = FALSE)

# Regions in a west-to-east-ish order, and populations within a region ordered
# by longitude, so that the population axis of every plot follows geography.
region_order <- c("AFRICA", "MIDDLE_EAST", "EUROPE", "CENTRAL_SOUTH_ASIA",
                  "EAST_ASIA", "OCEANIA", "AMERICA")
stopifnot(setequal(region_order, unique(m$region)))

pop_tab <- unique(m[, c("population", "region", "longitude", "latitude")])
pop_tab <- pop_tab[order(match(pop_tab$region, region_order), pop_tab$longitude,
                         pop_tab$population), ]
pop_order <- pop_tab$population

meta <- data.frame(
  sample    = m$sample,
  pop       = factor(m$population, levels = pop_order),
  region    = factor(m$region,     levels = region_order),
  latitude  = m$latitude,
  longitude = m$longitude,
  sex       = m$sex,
  stringsAsFactors = FALSE)
stopifnot(!any(is.na(meta$pop)), !any(is.na(meta$region)))

# ---- haplotype matrix -------------------------------------------------------

samples <- readLines(paste0(prefix, ".samples"))
snps    <- read.delim(paste0(prefix, ".snps"), header = FALSE,
                      col.names = c("snp", "ref", "alt", "aa", "daf"),
                      stringsAsFactors = FALSE)
ns <- length(samples)
nv <- nrow(snps)

# The .hap file is one SNP per line: 2*ns characters '0'/'1' then a newline, so
# reading it as raw bytes and reshaping gives a (2*ns+1) x nv array whose last
# row is the newlines. Dropping that row leaves haplotypes in rows and SNPs in
# columns, which is the orientation the fit wants. Stored as raw (1 byte per
# entry) to keep the object near 300 Mb rather than 2.5 Gb.
hapf <- paste0(prefix, ".hap")
w    <- 2L * ns + 1L
stopifnot(file.size(hapf) == as.double(w) * nv)
b <- readBin(hapf, "raw", n = file.size(hapf))
dim(b) <- c(w, nv)
stopifnot(all(b[w, 1:100] == as.raw(10)))     # newline in the last row
Y <- b[-w, , drop = FALSE]
rm(b); gc()

# Haplotype ids: sample_1 / sample_2, in the VCF's sample order.
rownames(Y) <- paste0(rep(samples, each = 2), "_", 1:2)
colnames(Y) <- snps$snp

stopifnot(all(samples %in% meta$sample))
meta <- meta[match(samples, meta$sample), ]
rownames(meta) <- NULL

saveRDS(meta, "data/hgdp_meta.rds")
saveRDS(list(Y = Y, snps = snps, samples = samples),
        "data/hgdp_hap_matrix.rds")

cat(sprintf("%d samples (%d haplotypes), %d SNPs\n", ns, 2 * ns, nv))
