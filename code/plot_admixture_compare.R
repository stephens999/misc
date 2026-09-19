# Structure plots for the two ADMIXTURE runs and the fastICA memberships they
# are being compared with, on a common component order and palette.
#
# Usage: Rscript code/plot_admixture_compare.R <outfile.png>

suppressMessages({library(ggplot2); library(cowplot); library(fastTopics)})
out <- commandArgs(trailingOnly = TRUE)[1]

cmp  <- readRDS("output/hgdp_admixture_compare.rds")
meta <- readRDS("data/hgdp_meta.rds")
pop  <- factor(as.character(cmp$pop), levels = levels(meta$pop))

sets <- list(`fastICA memberships (K = 25 maxima)` = cmp$Q_fastica,
             `ADMIXTURE, random start`             = cmp$Q_random,
             `ADMIXTURE, fastICA-initialized`      = cmp$Q_init)
cols <- Polychrome::glasbey.colors(ncol(cmp$Q_fastica) + 1)[-1]

mk <- function(nm, showx) {
  Q <- sets[[nm]]
  colnames(Q) <- make.unique(cmp$labels)
  set.seed(1)
  p <- structure_plot(Q, grouping = pop, gap = 8, verbose = FALSE,
                      colors = cols) +
    labs(y = "membership", fill = NULL, title = nm) +
    theme(legend.position = "none",
          plot.title = element_text(size = 9))
  if (showx) p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                                  hjust = 1, size = 5.5))
  else p + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
}
ps <- Map(mk, names(sets), c(FALSE, FALSE, TRUE))
leg <- get_legend(
  structure_plot(`colnames<-`(cmp$Q_fastica, make.unique(cmp$labels)),
                 grouping = pop, gap = 8, verbose = FALSE, colors = cols) +
    labs(fill = NULL) +
    theme(legend.position = "bottom", legend.text = element_text(size = 5.5),
          legend.key.size = unit(0.28, "cm")) +
    guides(fill = guide_legend(nrow = 3)))

g <- plot_grid(plot_grid(plotlist = ps, ncol = 1, align = "v",
                         rel_heights = c(1, 1, 1.5)),
               leg, ncol = 1, rel_heights = c(1, 0.16))
ggsave(out, g, width = 12, height = 10.5, dpi = 115)
cat("wrote", out, "\n")
