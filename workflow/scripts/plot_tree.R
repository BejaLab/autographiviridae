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
library(ggstar)
library(castor)

with(snakemake@input, {
    tree_file <<- tree
    ref_file <<- xlsx
    imgvr_file <<- metadata
})
with(snakemake@params, {
    ingroup_group <<- ingroup
    outgroup_group <<- outgroup
})

tree_file <- "analysis/phylophlan_markers/phylophlan.tre.treefile"
ref_file <- "input/phages.xlsx"
imgvr_file <- "analysis/exonuclease/IMGVR_metadata.csv"
cluster_file <- "input/host_clusters.xlsx"
strain_file <- "input/host_strains.xlsx"
exo_nbla_files <- Sys.glob("analysis/exonuclease/*_nblA.tab")
single_gene_blast_files <- Sys.glob("analysis/contigs/*_isolated_genes.outfmt6")

plot_file <- "tmp.svg"

ingroup_group <- "cyanophage"
outgroup_group <- c("pelagiphage", "coliphage")

strains <- read_excel(strain_file) %>%
   rename(Host = Strain, Host.taxonomy.prediction = Taxonomy)
clusters <- read_excel(cluster_file) %>%
   rename(Host.genus = Genus, Host.cluster = Cluster, Host.color = color)

col_names <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "stitle")
single_gene_blast <- lapply(single_gene_blast_files, read.table, col.names = col_names, sep = "\t") %>%
   bind_rows %>%
   arrange(-bitscore) %>%
   extract(sseqid, into = c("Phage", "gene"), regex = "(.+)_([^_]+)") %>%
   extract(stitle, into = "Host", regex = "[^ ]+ ([^,]+)") %>%
   distinct(Phage, .keep_all = T) %>%
   filter(pident > 95) %>%
   select(Accession = qseqid, Phage, Host, pident) %>%
   mutate(Host.evidence = "blast") %>%
   left_join(strains, by = "Host")

col_names <- c("Accession", "Ref", "Ref_clade", "Ref_identity", "Ref_length", "nblA_ORF", "Score", "Evalue")
nbla <- lapply(exo_nbla_files, read.table, header = F, row.names = NULL, col.names = col_names, sep = "\t") %>%
   bind_rows %>%
   filter(Score > 0) %>%
   select(Accession) %>%
   mutate(nblA = "exonuclease")

ref_metadata <- read_excel(ref_file) %>%
    left_join(strains, by = "Host")
ingroup <- filter(ref_metadata, Group %in% ingroup_group) %>%
    pull(Accession)
outgroup <- filter(ref_metadata, Group %in% outgroup_group) %>%
    pull(Accession)

root_tree <- function(tree, some_outgroup, some_ingroup) {
    tree <- root(tree, some_outgroup[1], resolve.root = T)
    ingroup_node <- getMRCA(tree, some_ingroup)
    tree <- root(tree, node = ingroup_node, resolve.root = T)
    ingroup_node <- getMRCA(tree, some_ingroup)
    extract.clade(tree, ingroup_node)
}
get_taxon <- function(.data, rank, col_prefix = "Host.") {
    regex <- sprintf("%s__([^;]+);", substr(rank, 1, 1))
    col_name <- sprintf("%s%s", col_prefix, rank)
    extract(.data, Host.taxonomy.prediction, into = col_name, regex = regex, remove = F)
}
get_mrca <- function(phylo, tips) {
    getMRCA(phylo, tips) %>%
        replace(is.null(.), NA)
}
add_mrca <- function(tree, colname, max_dist = 0.8) {
    colname <- deparse(substitute(colname))
    phylo <- to_treedata(tree)@phylo
    tree_df <- mutate(tree, my_column = !!as.name(colname))
    lumped <- distTips(phylo) %>%
        as.matrix %>%
        as.data.frame %>%
        rownames_to_column("A") %>%
        gather(B, dist, -A) %>%
        filter(dist <= max_dist) %>%
        left_join(tree_df, by = c(A = "label")) %>%
        filter(!is.na(my_column)) %>%
        with(setNames(my_column, B))
    mrca <- mutate(tree, my_column = lumped[label]) %>%
        group_by(my_column) %>%
        mutate(is.tip = label %in% phylo$tip.label) %>%
        mutate(no_data = all(is.na(my_column))) %>%
        mutate(mrca = get_mrca(phylo, node[is.tip])) %>%
        mutate(mrca = ifelse(no_data | is.na(mrca), node, mrca)) %>%
        group_by(mrca) %>%
        mutate(enough_tips = sum(is.tip) > 1) %>%
        mutate(ifelse(node == mrca & enough_tips, first(na.omit(my_column)), NA)) %>%
        pull
    tree[[paste0(colname, "_mrca")]] <- mrca
    return(tree)
}
to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}
parse_support <- function(tree) {
    mutate(tree, support = ifelse(node %in% parent, label, NA)) %>%
        separate(support, into = c("SH_aLRT", "UFBoot"), sep = "/", convert = T)
}
pick_support_values <- function(tree, support_col, max_dist = 1) {
    support_col <- deparse(substitute(support_col))
    phylo <- as.treedata(tree)@phylo
    parents <- get_all_distances_to_root(phylo) %>%
        {data.frame(dist = ., node = 1:length(.))} %>%
        filter(dist <= max_dist, node > length(phylo$tip.label)) %>%
        pull(node)
    mutate(tree, support = ifelse(parent %in% parents, !!as.name(support_col), NA))
}
tree <- read.tree(tree_file) %>%
    root_tree(outgroup, ingroup)

imgvr_metadata <- read.csv(imgvr_file) %>%
    mutate(Accession = paste0(UVIG, "_1")) %>%
    mutate(Host.unknown = Host.prediction.method == '' | is.na(Host.prediction.method)) %>%
    left_join(single_gene_blast, by = "Accession") %>%
    mutate(Host.taxonomy.prediction = ifelse(!Host.unknown & is.na(Host.taxonomy.prediction.y), Host.taxonomy.prediction.x, Host.taxonomy.prediction.y)) %>%
    mutate(Host.evidence = case_when(Host.unknown ~ Host.evidence, Topology == "Provirus" ~ "Prophage", T ~ Host.prediction.method))

gov2_metadata <- data.frame(Accession = tree$tip.label) %>%
    filter(grepl("^Station", Accession)) %>%
    mutate(Ecosystem.classification = "Environmental;Aquatic;Marine;") %>%
    left_join(single_gene_blast, by = "Accession")

env_metadata <- bind_rows(imgvr_metadata, gov2_metadata)

metadata <- bind_rows(ref_metadata, env_metadata) %>%
    distinct(Accession, .keep_all = T) %>%
    mutate(Ecosystem = sub(";[^;]*$", "", Ecosystem.classification)) %>%
    mutate(Ecosystem = ifelse(Ecosystem == ";;", NA, Ecosystem)) %>%
    mutate(Prophage = Topology == "Provirus") %>%
    get_taxon("phylum") %>%
    get_taxon("family") %>%
    get_taxon("genus") %>%
    left_join(clusters, by = "Host.genus") %>%
    left_join(nbla, by = "Accession") %>%
    mutate(nblA = ifelse(is.na(nblA.x) | nblA.x == '', nblA.y, nblA.x)) %>%
    mutate(nblA = ifelse(is.na(nblA) | nblA == '', NA, nblA)) %>%
    select(Accession, Phage, Group, Clade, Subclade, Host, Host.phylum, Host.genus, Host.cluster, Prophage, Host.evidence, nblA, Ecosystem)

tree_tib <- as_tibble(tree) %>%
    parse_support %>%
    pick_support_values(UFBoot, max_dist = 5) %>%
    left_join(metadata, by = c(label = "Accession")) %>%
    add_mrca(Clade) %>%
    add_mrca(Subclade, max_dist = 0.745)

evidence.shapes <- c(
    `Isolate taxonomy` = 15,
    `Prophage` = 0,
    `blast` = 2,
    `CRISPR spacer match` = 1,
    `k-mer match` = 1
)
cluster.colors <- with(clusters, setNames(Host.color, Host.cluster))
ecosystems <- select(metadata, Accession, Ecosystem_tile = Ecosystem)
ecosystem.colors <- c(
    `Environmental;Aquatic;Marine` = "blue",
    `Environmental;Aquatic;Freshwater` = "cyan"
)

p <- ggtree(to_treedata(tree_tib), size = 0.1, layout = "rectangular", open.angle = 45) +
    geom_nodepoint(aes(x = branch, subset = !is.na(support) & support >= 95), size = 0.2, color = "#4d4d4dff") +
    geom_tiplab(aes(subset = !is.na(Phage), label = Phage), size = 2, align = T, linesize = 0.1) +
    geom_tippoint(aes(subset = !is.na(nblA), color = nblA), size = 0.5, color = "red") + new_scale_color() +
    geom_tippoint(aes(subset = !is.na(Host.cluster), x = x + 0.05, color = Host.cluster, shape = Host.evidence)) +
    scale_color_manual(values = cluster.colors) +
    scale_shape_manual(values = evidence.shapes) +
    geom_fruit(aes(x = 0, y = Accession, fill = Ecosystem_tile), ecosystems, geom_tile, axis.params = c(axis = "x", text = "Ecosystem"), offset = 0.5) +
    scale_fill_manual(values = ecosystem.colors) + new_scale_fill() +
    geom_cladelab(mapping = aes(subset = !is.na(Subclade_mrca), node = node, label = Subclade_mrca), align = T, offset = 0.3) +
    geom_treescale(width = 0.4)

ggsave(plot_file, p, width = 8, height = 12)
