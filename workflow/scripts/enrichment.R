library <- function(...) suppressPackageStartupMessages(base::library(...))
library(ggplot2)
library(tibble)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(enrichplot)
library(stringr)
library(shadowtext)
library(readxl)

if (interactive()) {
    results_file <- "input/Day2_NblAOE_For_Go Enrichment.xlsx"
    annot_file <- "analysis/enrichment/go-terms-MF.tsv"
    output_file <- "tmp.pdf"
    sheet <- "SemiN Ratio <-1 only NblA"
    alpha <- 0.05
} else {
    with(snakemake@input, {
        results_file <<- results
        annot_file <<- annotation
    })
    output_file <- unlist(snakemake@output)
    with(snakemake@params, {
        alpha <<- alpha
        sheet <<- sheet
    })
}
ora_generic <- function(gene_set, term2name) {
    gene_universe <- unique(term2gene$gene)
    term2name <- distinct(term2name, annot, description)
    enricher(gene_set, universe = gene_universe, pAdjustMethod = "BH", qvalueCutoff = alpha, TERM2GENE = term2gene, TERM2NAME = term2name)
}

term2gene <- read.table(annot_file, col.names = c("annot", "gene", "description"), sep = "\t", quote = "", fill = T) %>%
    mutate(description = ifelse(is.na(description), annot, description))
data <- read_excel(results_file, sheet = sheet)
gene_set <- filter(data, ACC %in% term2gene$gene) %>%
    pull(ACC)

enrichment <- ora_generic(gene_set, term2gene)
results <- enrichment@result %>%
    filter(p.adjust <= alpha, Count > 1) %>%
    separate_rows(geneID) %>%
    mutate(n_reg = n()) %>%
    ungroup %>%
    separate(GeneRatio, into = c("GeneCount", "GeneTotal"), convert = T) %>%
    mutate(GeneRatio = GeneCount / GeneTotal)

p <- ggplot(results, aes(x = GeneRatio, y = Description, color = p.adjust, size = Count)) +
    geom_point() +
    geom_shadowtext(aes(label = Count), bg.r = 0.02, color = "white") +
    scale_color_gradient(name = "p adjusted", limits = c(0, alpha), low = "red", high = "blue") +
    scale_x_continuous(name = "Gene ratio", label = abs, expand = expansion(mult = 0.2)) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) +
    scale_size(range = c(3, 10)) +
    guides(size = "none") +
    theme_bw() +
    theme(axis.title.y = element_blank())
num_cats <- nrow(distinct(all_results, ID))
num_stages <- nrow(distinct(all_results, stage))
ggsave(output_file, p, width = num_stages * 1.8 + 2, height = num_cats * 0.3 + 1)
