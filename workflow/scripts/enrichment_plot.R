library <- function(...) suppressPackageStartupMessages(base::library(...))
library(ggplot2)
library(gridExtra)
library(grid)
library(dplyr)
library(tidyr)
library(shadowtext)
library(stringr)

csv_files <- unlist(snakemake@input)
output_file <- unlist(snakemake@output)
alpha <- unlist(snakemake@params["alpha"])

results <- lapply(csv_files, read.csv) %>%
    bind_rows %>%
    mutate(EnrichmentFold = GeneRatio / BgRatio) %>%
    mutate(Description = str_replace_all(Description, " ", "\\u00a0"))

p <- ggplot(results, aes(x = EnrichmentFold, y = Description, color = p.adjust, size = Count)) +
    facet_grid(ONTOLOGY ~ day, scales = "free_y", space = "free_y") +
    geom_point() +
    geom_shadowtext(aes(label = Count), bg.r = 0.02, color = "white") +
    scale_color_gradient(name = "p adjusted", limits = c(0, alpha), low = "red", high = "blue") +
    scale_x_continuous(name = "Enrichment fold", label = abs, expand = expansion(mult = 0.2)) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 30)) +
    scale_size(range = c(3, 10)) +
    guides(size = "none") +
    theme_bw() +
    theme(axis.title.y = element_blank(), strip.background = element_blank(), strip.text.y.right = element_text(angle = 0), axis.text.y = element_text(size = 8)) +
    labs(color = "adjusted p")
height <- 1 + 0.25 * nrow(results)
ggsave(output_file, p, width = 15, height = height)
