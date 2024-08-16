library <- function(...) suppressPackageStartupMessages(base::library(...))
library(AnnotationForge)
library(AnnotationDbi)
library(XML)
library(dplyr)
library(tidyr)

uniprot_file_path <- unlist(snakemake@input)
orgdb_dir <- unlist(snakemake@output)
# Parsing the XML for Uniprot IDs and GO accessions
uniprot_xml <- xmlTreeParse(uniprot_file_path, useInternalNodes = T)
entries <- getNodeSet(uniprot_xml, "//ns:entry", namespace = c(ns = "http://uniprot.org/uniprot"))

uniprot_ids <- sapply(entries, function(entry) {
  return(xmlValue(entry[["accession"]]))
})
uniprot_genes <- sapply(entries, function(entry) {
  return(xmlValue(entry[["gene"]][["name"]][[1]]))
})
uniprot_names <- sapply(entries, function(entry) {
  return(xmlValue(entry[["protein"]][["recommendedName"]][["fullName"]]))
})
gene_info <- data.frame(GID = uniprot_ids, SYMBOL = uniprot_genes, GENENAME = uniprot_names) %>%
  replace_na(list(GENENAME = ""))

go_terms <- lapply(entries, function(entry) {
  go_nodes <- getNodeSet(entry, ".//ns:dbReference[@type='GO']", namespace = c(ns = "http://uniprot.org/uniprot"))
  sapply(go_nodes, function(go_node) {
    return(xmlGetAttr(go_node, "id"))
  })
})
go_df <- setNames(go_terms, uniprot_ids) %>%
  lapply(data.frame) %>%
  bind_rows(.id = "GID") %>%
  setNames(c("GID", "GO")) %>%
  mutate(EVIDENCE = "IEA")
chr <- data.frame(GID = uniprot_ids, CHROMOSOME = 1)

organism <- entries[[1]][["organism"]]
binomen <- strsplit(xmlValue(organism[["name"]]), " ")
genus <- binomen[[1]][1]
species <- paste(binomen[[1]][3:4], collapse = "")
tax_id <- xmlGetAttr(organism[["dbReference"]], "id")
author <- "Author <example@server.com>"

makeOrgPackage(
  gene_info = gene_info,
  go = go_df,
  goTable = "go",
  chromosome = chr,
  version = "0.1",
  author = author, maintainer = author,
  outputDir = dirname(orgdb_dir),
  tax_id = tax_id,
  genus = genus,
  species = species
)
# Sys.chmod(orgdb_file, "666", use_umask = F)
