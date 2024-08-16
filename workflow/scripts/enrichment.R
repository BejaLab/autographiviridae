library <- function(...) suppressPackageStartupMessages(base::library(...))
library(tibble)
library(dplyr)
library(tidyr)
library(clusterProfiler)
library(readxl)
library(AnnotationDbi)

if (interactive()) {
    results_file <- "input/Day2_NblAOE_For_Go Enrichment.xlsx"
    output_file <- "tmp.csv"
    sheet <- "SemiN Ratio <-1 only NblA"
    orgdb_file <- "analysis/orgDb/org.SWH8109.eg.db/inst/extdata/org.SWH8109.eg.sqlite"
    alpha_threshold <- 0.05
    cutoff <- 0.5
} else {
    with(snakemake@input, {
        results_file <<- results
        orgdb_dir <<- orgdb
    })
    orgdb_file <- Sys.glob(file.path(orgdb_dir, "inst", "extdata", "*.sqlite"))
    orgdb_package <- basename(orgdb_dir)
    output_file <- unlist(snakemake@output)
    with(snakemake@params, {
        alpha_threshold <<- alpha
        cutoff <<- cutoff
    })
    with(snakemake@wildcards, {
        day <<- day
        sheet <<- sheet
    })
}
orgdb <- loadDb(orgdb_file)

data <- read_excel(results_file, sheet = sheet) %>%
    dplyr::rename(protein = 1)

enrichment <- enrichGO(data$protein, orgdb, keyType = "GID", ont = "ALL", pAdjustMethod = "BH", qvalueCutoff = alpha_threshold)
sim_enrich <- simplify(enrichment, cutoff = cutoff)

results <- mutate(sim_enrich@result, sheet = sheet, day = day) %>%
    filter(p.adjust <= alpha_threshold) %>%
    group_by(ONTOLOGY, geneID) %>%
    arrange(p.adjust) %>%
    mutate(Description = paste(Description, collapse = '/')) %>%
    distinct(ONTOLOGY, geneID, .keep_all = T) %>%
    ungroup %>%
    separate(GeneRatio, into = c("GeneCount", "GeneTotal"), convert = T) %>%
    separate(BgRatio, into = c("BgCount", "BgTotal"), convert = T) %>%
    mutate(GeneRatio = GeneCount / GeneTotal, BgRatio = BgCount / BgTotal)
write.csv(results, file = output_file, row.names = F)
