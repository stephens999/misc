# Nearest-gene annotation for the TGP SNPs, for labelling outliers in
# analysis/fastica_1kg_structure.Rmd.
#
# The 1000 Genomes phase 3 coordinates are GRCh37, so we use the UCSC hg19
# refGene table. (The installed org.Hs.eg.db is GRCh38: it puts LCT at
# 135.79 Mb where GRCh37 has it at 136.55 Mb, so it cannot be used here.)
# The table is downloaded on first run and cached in data/, which is gitignored;
# the small derived annotation is what gets committed.
#
# Usage: Rscript code/annotate_snps.R
# Writes output/tgp_snp_genes.rds: one row per SNP with its nearest gene and
# the distance to that gene's transcript bounds (0 if inside).

ref_gz <- "data/refGene_hg19.txt.gz"
if (!file.exists(ref_gz))
  download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refGene.txt.gz",
                ref_gz, quiet = TRUE)

# refGene columns: bin, name, chrom, strand, txStart, txEnd, ... , name2 (symbol)
rg <- read.delim(gzfile(ref_gz), header = FALSE, stringsAsFactors = FALSE)[, c(3, 5, 6, 13)]
names(rg) <- c("chrom", "start", "end", "gene")
rg <- rg[grepl("^chr([0-9]+|X|Y)$", rg$chrom), ]
rg$chrom <- sub("^chr", "", rg$chrom)

# One interval per gene per chromosome, spanning all its transcripts.
rg <- do.call(rbind, lapply(split(rg, list(rg$chrom, rg$gene), drop = TRUE),
                            function(d) data.frame(chrom = d$chrom[1],
                                                   start = min(d$start),
                                                   end   = max(d$end),
                                                   gene  = d$gene[1])))
cat("genes:", nrow(rg), "on", length(unique(rg$chrom)), "chromosomes\n")

snps <- colnames(readRDS("output/tgp_B_pruned_maxima.rds"))
chr  <- sub(":.*", "", snps)
pos  <- as.numeric(sub(".*:", "", snps))

# Nearest gene within each chromosome: distance is 0 inside a gene, otherwise
# the gap to the closest transcript bound.
out <- data.frame(snp = snps, gene = NA_character_, dist = NA_real_,
                  stringsAsFactors = FALSE)
for (cc in unique(chr)) {
  g <- rg[rg$chrom == cc, ]
  i <- which(chr == cc)
  if (!nrow(g)) next
  d <- outer(pos[i], g$start, function(p, s) pmax(s - p, 0)) +
       outer(pos[i], g$end,   function(p, e) pmax(p - e, 0))
  k <- apply(d, 1, which.min)
  out$gene[i] <- g$gene[k]
  out$dist[i] <- d[cbind(seq_along(i), k)]
}
cat("annotated:", sum(!is.na(out$gene)), "of", nrow(out), "SNPs |",
    "within a gene:", sum(out$dist == 0, na.rm = TRUE), "\n")

saveRDS(out, "output/tgp_snp_genes.rds")
