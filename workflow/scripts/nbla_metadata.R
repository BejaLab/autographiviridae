library <- function(...) suppressPackageStartupMessages(base::library(...))
library(dplyr)
library(tidyr)
library(treeio)

with(snakemake@input, {
    jtree_file <<- jtree
    ingroup_file <<- ingroup
    pro_syn_file <<- pro_syn
    cyanorak_file <<- cyanorak
    img_file <<- img
    refs_file <<- refs
    cyanobacteria_file <<- cyanobacteria
})
output_file <- unlist(snakemake@output)

with(snakemake@params, {
    show_genomes <<- show_genomes
})

read_picocyano <- function(file_name) {
    synonyms <- c(LLII = "LLII/III", LLIII = "LLII/III")
    read.csv(file_name) %>%
        select(genome = Name, SubCluster, Clade, SubClade) %>%
        mutate(ProSynCya = sub("_.+", "", genome)) %>%
        mutate(Clade = recode(Clade, !!!synonyms), SubClade = recode(SubClade, !!!synonyms)) %>%
        mutate(clade = case_when(
            ProSynCya == "Pro" ~ "Pro",
            ProSynCya == "Syn" | ProSynCya == "Cya" ~ "Syn"
        )) %>%
        mutate(picocyano_group = case_when(
            clade == "Pro" ~ paste("Prochlorococcus subclade", SubClade),
            clade == "Syn" ~ paste("\"Synechococcus\" subcluster", SubCluster)
        )) %>%
        mutate(alias = case_when(clade == "Pro" ~ genome))
}

read_fasta <- function(filename) {
    read.fasta(filename) %>%
        {data.frame(desc = names(.))} %>%
        mutate(label = sub(" .+", "", desc))
}
genomes <- read.jtree(jtree_file) %>%
    as_tibble %>%
    filter(!is.na(Accession), !is.na(nblA)) %>%
    mutate(nbla_category = paste0("T7_", Clade_assigned)) %>%
    mutate(nbla_category2 = paste0("T7_", nblA)) %>%
    select(genome = Accession, nbla_category2, nbla_category, alias = Phage)

ingroup <- read_fasta(ingroup_file) %>%
    extract(desc, into = "Accession", regex = "cds:([^: ]+)") %>%
    mutate(genome = sub("_[0-9]+$", "", label)) %>%
    left_join(genomes, by = "genome") %>%
    mutate(taxon = "Cyanophages")

picocyano <- bind_rows(
    read_picocyano(cyanorak_file),
    read_picocyano(img_file)
)
pro_syn <- read_fasta(pro_syn_file) %>%
    mutate(genome = sub("^([^_]+_[^_]+(_[A-Z][0-9]+)?)_.+", "\\1", label)) %>%
    left_join(picocyano, by = "genome") %>%
    mutate(alias = ifelse(genome %in% show_genomes, sub("^[A-Za-z]+_", "", genome), NA)) %>%
    mutate(taxon = "Cyanobacteria", nbla_category2 = clade, nbla_category = picocyano_group)

refs <- read_fasta(refs_file) %>%
    extract(desc, into = "taxon", regex = "Tax=([^=]+)(?=\\s\\w+=|$)", remove = F) %>%
    extract(desc, into = "alias", regex = "Alias=([^=]+)(?=\\s\\w+=|$)", remove = F) %>%
    extract(desc, into = "nbla_category",  regex = "Category=([^=]+)(?=\\s\\w+=|$)", remove = F) %>%
    extract(desc, into = "nbla_category2", regex = "Category2=([^=]+)(?=\\s\\w+=|$)", remove = F)

cyanobacteria <- read_fasta(cyanobacteria_file) %>%
    mutate(taxon = "Cyanobacteria", nbla_category = "Other_cyanobacteria")

metadata <- bind_rows(list(ingroup, pro_syn, refs, cyanobacteria)) %>%
    select(label, taxon, nbla_category, nbla_category2, alias) %>%
    mutate(nbla_category2 = ifelse(is.na(nbla_category2), nbla_category, nbla_category2))

write.csv(metadata, file = output_file, quote = F, row.names = F, na = "")
