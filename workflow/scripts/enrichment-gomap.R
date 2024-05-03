library <- function(...) suppressPackageStartupMessages(base::library(...))
library(dplyr)
library(ontologyIndex)
library(clusterProfiler)

with(snakemake@input, {
    tsv_file <<- tsv
    obo_file <<- obo
})
output_file <- snakemake@output
obo <- get_ontology(obo_file, extract_tags = "everything")
term2gene <- read.table(tsv_file, col.names = c("go_id", "gene"))

ontologies <- list(
    BP = "biological_process",
    CC = "cellular_component",
    MF = "molecular_function"
)
descriptions <- data.frame(GO = names(obo$name), description = unname(obo$name))

gomap <- mutate(term2gene, ont = unlist(obo$namespace[go_id])) %>%
    split(f = .$ont) %>%
    lapply(select, annot = 1, gene = 2) %>%
    lapply(buildGOmap) %>%
    lapply(left_join, descriptions, by = "GO")

lapply(names(ontologies), function(ont) {
    full <- ontologies[[ont]]
    write.table(gomap[[full]], file = output_file[[ont]], sep = "\t", quote = F, col.names = F, row.names = F)
}) %>% invisible()
