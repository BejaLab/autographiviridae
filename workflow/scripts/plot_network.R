library <- function(...) suppressPackageStartupMessages(base::library(...))

library(treeio)
library(ggtree)
library(dplyr)
library(tidyr)
library(readxl)
library(phangorn)
library(ggnewscale)
library(ggplot2)
library(tanggle)
library(stringr)

with(snakemake@input, {
    network_file <<- network   # analysis/nbla/nbla_all_cdhit.nex
    cyanorak_file <<- cyanorak # analysis/cyanorak/organisms.csv
    cyanobact_file <<- cyanobacteria # analysis/nbla/nblA_cyanobacteria.faa
    phage_jtree_file <<- phage_jtree # output/phylophlan_markers9.jtree
    clstr_file <<- clstr # analysis/nbla/nblA_all_cdhit.faa.clstr
    refs_file <<- refs # input/profile/nblA_refs.faa
})
plot_file <- unlist(snakemake@output)

read_clstr <- function(file_name) {
    read.table(file_name, col.names = c("col1", "col2", "col3", "col4", "col5"), fill = T) %>%
        mutate(is.representative = col4 == "*") %>%
        mutate(cluster = ifelse(col1 == ">Cluster", col2, NA)) %>%
        fill(cluster) %>%
        mutate(label = str_sub(col3, 2, -4)) %>%
        filter(label != "") %>%
        group_by(cluster) %>%
        mutate(representative = label[is.representative]) %>%
        select(label, representative, cluster)
}
read_cyanorak <- function(file_name) {
    read.csv(cyanorak_file) %>%
        select(Prefix = Name, SubCluster, Clade, SubClade, Pigment.type) %>%
        mutate(ProSynCya = sub("_.+", "", Prefix)) %>%
        mutate(cyanorak_group = case_when(
            ProSynCya == "Pro" ~ paste("Prochlorococcus subclade", SubClade),
            ProSynCya == "Syn" | ProSynCya == "Cya" ~ paste("\"Synechococcus\" subcluster", SubCluster)
        ))
}
read_phages <- function(jtree_file) {
    as_tibble(read.jtree(phage_jtree_file)) %>%
        select(Accession, nblA, Phage, Host.cluster, Clade_assigned) %>%
        mutate(phages_group = paste("Autographiviridae clade", Clade_assigned))
}
read_fasta <- function(fasta_file) {
    read.fasta(fasta_file) %>%
        {data.frame(desc = names(.))} %>%
        separate(desc, into = c("label", "def"), sep = " ", extra = "merge") %>%
        extract(def, into = "Tax", regex = "Tax=([^=]+)(?=\\s\\w+=|$)", remove = F) %>%
        extract(def, into = "Alias", regex = "Alias=([^=]+)(?=\\s\\w+=|$)", remove = F)
}

network <- read.nexus.networx(network_file)
refs <- read_fasta(refs_file)
cyanorak <- read_cyanorak(cyanorak_file)
cyanobact <- read_fasta(cyanobact_file)
phages <- read_phages(phage_jtree_file)

clstr <- read_clstr(clstr_file)

node_data <- read_clstr(clstr_file) %>%
    mutate(Accession = sub("_[0-9]+$", "", label)) %>%
    mutate(Prefix = sub("^([^_]+_[^_]+).+", "\\1", Accession)) %>%
    left_join(cyanorak, by = "Prefix") %>%
    left_join(phages, by = "Accession") %>%
    left_join(refs, by = "label") %>%
    mutate(alias = case_when(
        !is.na(Alias) ~ Alias,
        !is.na(Phage) ~ Phage
    )) %>%
    mutate(Category = case_when(
        !is.na(phages_group) ~ "Cyanophages",
        !is.na(cyanorak_group) ~ "Cyanobacteria",
        Tax == "Cyanophages" ~ "Cyanophages",
        Tax == "Cyanobacteriota" ~ "Cyanobacteria",
        label %in% cyanobact$label ~ "Cyanobacteria"
    )) %>%
    mutate(Group = case_when(
        !is.na(phages_group) ~ phages_group,
        !is.na(cyanorak_group) ~ cyanorak_group,
        Tax == "Cyanophages" ~ "Other cyanophages",
        Tax == "Cyanobacteriota" ~ "Other cyanobacteria",
        label %in% cyanobact$label ~ "Other cyanobacteria"
    )) %>%
    group_by(label = representative) %>%
    summarize(alias = paste(na.omit(alias), collapse = "; "), Group = first(na.omit(Group)), Category = first(na.omit(Category)), .groups = "drop")

x_ext <- Vectorize(function(x0, y0, x1, y1, d) {
    dx <- x1 - x0
    dy <- y1 - y0
    return(x1 + dx * d / sqrt(dx^2 + dy^2))
}, vectorize.args = c("x0", "y0", "x1", "y1"))
y_ext <- Vectorize(function(x0, y0, x1, y1, d) {
    dx <- x1 - x0
    dy <- y1 - y0
    return(y1 + dy * d / sqrt(dx^2 + dy^2))
}, vectorize.args = c("x0", "y0", "x1", "y1"))

add_node_data <- function(p, more_data) {
    p$data <- left_join(p$data, more_data, by = "label")
    return(p)
}

p <- ggsplitnet(network, size = 0.1) %>% add_node_data(node_data) +
    geom_point2(aes(subset = !is.na(Group), x = x, y = y, color = Group, shape = Category)) +
    geom_tiplab2(aes(label = alias, x = x_ext(xend, yend, x, y, 0.01), y = y_ext(xend, yend, x, y, 0.01)), color = "gray30") +
    new_scale_color() +
    geom_treescale(width = 0.1) +
    coord_fixed()
ggsave(plot_file, p, width = 10, height = 8)
