library(dplyr)
library(tidyr)
library(treeio)
library(readxl)
library(ape)
library(ggnewscale)
library(ggtree)
library(ggtreeExtra)
library(tibble)
library(adephylo)
library(ggplot2)
library(castor)
library(jsonlite)
library(phangorn)

with(snakemake@input, {
    tree_file <<- tree
    cluster_file <<- clusters
    jtree_file <<- jtree
    fasta_file <<- fasta
})
with(snakemake@output, {
    plot_file <<- plot
    jtree_out_file <<- jtree
})
# tree_file <- "analysis/genomes/all_genomes_nblA.treefile"
# jtree_file <- "output/phylophlan_markers9.jtree"
# plot_file <- "tmp.svg"

clusters <- read_excel(cluster_file) %>%
    rename(Host.genus = Genus, Host.cluster = Cluster, Host.color = color)

fasta <- read.fasta(fasta_file) %>%
    as.character %>%
    {data.frame(label = names(.), sequence = sapply(., paste, collapse = ""))} %>%
    mutate(sequence = toupper(sequence)) %>%
    separate(label, into = c("label", "description"), sep = " ", extra = "merge")

metadata <- read.jtree(jtree_file) %>%
    as_tibble %>%
    filter(!is.na(Accession), !is.na(nblA)) %>%
    select(Accession, nblA, Phage, Host.cluster, Clade_assigned)

to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}   
tree <- read.tree(tree_file) %>%
    as_tibble %>%
    mutate(is.tip = ! node %in% parent) %>%
    mutate(UFBoot = as.numeric(ifelse(is.tip, NA, label))) %>%
    mutate(Accession = ifelse(is.tip, sub("_[0-9]+$", "", label), NA)) %>%
    left_join(metadata, by = "Accession") %>%
    left_join(fasta, by = "label") %>%
    to_treedata
write.jtree(tree, file = jtree_out_file)

cluster.colors <- with(clusters, setNames(Host.color, Host.cluster))
p <- ggtree(tree, aes(color = Clade_assigned), size = 0.2, layout = "ape", open.angle = 45) +
    new_scale_color() +
    geom_point2(aes(subset = UFBoot == 100, x = branch.x, y = branch.y), size = 0.2, color = "black") +
    geom_tiplab2(aes(subset = !is.na(Phage), label = Phage, color = Host.cluster), size = 2, fontface = "bold") +
    scale_color_manual(values = cluster.colors) + new_scale_color() +
    geom_tippoint(aes(subset = !is.na(nblA), color = nblA), size = 0.5) + new_scale_color() +
    geom_treescale(width = 0.4)
ggsave(plot_file, p, width = 10, height = 6)
