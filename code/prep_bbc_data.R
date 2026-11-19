# Build the BBC news document-term matrix used by analysis/fastica_bbc.Rmd.
#
# Weine et al. (arXiv 2601.05845) analyse this corpus in
# log1pNMF/inst/paper_figures/fig_scripts/bbc_fig.R, starting from a csv
# (bbc_news_text_complexity_summarization.csv) that is not distributed with the
# package. That csv is a deduplicated copy of the original BBC news dataset of
# Greene and Cunningham (2006), which is public:
#
#   http://mlg.ucd.ie/files/datasets/bbc-fulltext.zip
#
# The zip has 2225 articles in five category directories, of which 2127 are
# distinct as plain text. Dropping the 98 exact duplicates and running the
# tokenisation from bbc_fig.R verbatim reproduces their matrix exactly:
# 2127 x 5861, with every one of their 5861 terms matched.
#
# Usage: Rscript code/prep_bbc_data.R [zipfile]
# Writes data/bbc_dtm.rds: counts (2127 x 5861 dgCMatrix), labels, file names,
#        headlines.

suppressPackageStartupMessages({
  library(tm); library(SnowballC); library(Matrix)
})

# The only non-ASCII content in the corpus is the pound sign, as valid UTF-8
# (bytes c2 a3), in 647 of the 2225 files. In a non-UTF-8 locale -- which is
# what Rscript gets when LANG is unset -- tm escapes those bytes to the
# literal eight-character string "<c2><a3>", which then splits "£bn" into two
# junk terms and leaves the escape text in every table downstream. Forcing a
# UTF-8 locale here keeps "£bn" as one term and makes the vocabulary match
# Weine et al.'s exactly.
if (!isTRUE(l10n_info()$`UTF-8`)) {
  for (loc in c("en_US.UTF-8", "C.UTF-8", "UTF-8"))
    if (suppressWarnings(Sys.setlocale("LC_CTYPE", loc)) != "") break
  if (!isTRUE(l10n_info()$`UTF-8`))
    warning("could not set a UTF-8 locale; the pound sign will be mangled")
}

args <- commandArgs(trailingOnly = TRUE)
zipf <- if (length(args) > 0) args[1] else "data/bbc-fulltext.zip"
if (!file.exists(zipf)) {
  cat("downloading", zipf, "\n")
  download.file("http://mlg.ucd.ie/files/datasets/bbc-fulltext.zip", zipf,
                quiet = TRUE)
}
root <- file.path(tempdir(), "bbc-fulltext")
unzip(zipf, exdir = root)
root <- file.path(root, "bbc")

cats  <- c("business", "entertainment", "politics", "sport", "tech")
files <- unlist(lapply(cats, function(cc)
  file.path(cc, list.files(file.path(root, cc)))))
labs  <- factor(dirname(files), levels = cats,
                labels = c("Business", "Entertainment", "Politics", "Sports",
                           "Tech"))

# The first line of each file is the headline; bbc_fig.R's csv has headline and
# body in one `text` field, so paste the lines together.
txt <- vapply(file.path(root, files),
              function(f) {
                x <- paste(readLines(f, warn = FALSE), collapse = " ")
                Encoding(x) <- "UTF-8"
                x
              },
              character(1))
keep <- !duplicated(txt)
cat(sprintf("%d articles, %d distinct\n", length(txt), sum(keep)))
txt <- txt[keep]; files <- files[keep]; labs <- labs[keep]

# First line of each file is the headline. Kept so that analyses can name the
# documents a factor loads on, not just count them.
headlines <- vapply(file.path(root, files),
                    function(f) trimws(readLines(f, warn = FALSE)[1]),
                    character(1), USE.NAMES = FALSE)

# ---- tokenisation, exactly as in bbc_fig.R ---------------------------------

cp <- VCorpus(VectorSource(txt))
addspace <- content_transformer(function(x, pattern) gsub(pattern, " ", x))
cp <- tm_map(cp, addspace, "-")
cp <- tm_map(cp, removeNumbers)
cp <- tm_map(cp, content_transformer(tolower))
cp <- tm_map(cp, removeWords, stopwords("SMART"))
cp <- tm_map(cp, removePunctuation)
cp <- tm_map(cp, stripWhitespace)
cp <- tm_map(cp, stemDocument)

dtm <- DocumentTermMatrix(cp)
X   <- sparseMatrix(i = dtm$i, j = dtm$j, x = dtm$v,
                    dims = c(dtm$nrow, dtm$ncol))
colnames(X) <- dtm$dimnames$Terms
rownames(X) <- files

# Keep terms appearing in more than four documents, as they do.
X <- X[, Matrix::colSums(X > 0) > 4, drop = FALSE]
cat(sprintf("%d documents x %d terms, %.2f%% non-zero\n",
            nrow(X), ncol(X), 100 * length(X@x) / prod(dim(X))))
print(table(labs))

saveRDS(list(X = X, labels = labs, files = files, headlines = headlines,
             size = Matrix::rowSums(X)),
        "data/bbc_dtm.rds")
cat("wrote data/bbc_dtm.rds\n")
