library <- function(...) suppressPackageStartupMessages(base::library(...))

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
    ref_file <<- refs
    imgvr_file <<- imgvr
    cluster_file <<- clusters
    strain_file <<- strains
    exo_nbla_files <<- exo_nbla
    psba_file <<- psba
    pol_file <<- pol
    genome_nbla_file <<- genome_nbla
})
with(snakemake@output, {
    plot_all_file <<- plot_all
    plot_ref_file <<- plot_ref
    jtree_file <<- jtree
})
with(snakemake@params, {
    w_all <<- w_all
    h_all <<- h_all
    w_ref <<- w_ref
    h_ref <<- h_ref
    nbla_min_score <<- nbla_min_score
    max_mrca_dist <<- max_mrca_dist
    pol_domain <<- pol_domain
    exo_domain <<- exo_domain
})
# tree_file <- "analysis/phylophlan_markers/phylophlan.tre.treefile"
# ref_file <- "input/phages.xlsx"
# imgvr_file <- "analysis/exonuclease/IMGVR_metadata.csv"
# cluster_file <- "input/host_clusters.xlsx"
# strain_file <- "input/host_strains.xlsx"
# exo_nbla_files <- Sys.glob("analysis/exonuclease/*_nblA.tab")
# genome_nbla_file <- "analysis/genomes/nblA.jsonl"
# plot_file <- "tmp.svg"
# jtree_file <- "tmp.jtree"

read_jsonl <- function(filename) {
    readLines(filename) %>%
        lapply(fromJSON) %>%
        lapply(unlist) %>%
        bind_rows
}
read_tblout <- function(filename) {
    col_names <- c("target_name", "target_accession", "query_name", "query_accession", "full_E_value", "full_score", "full_bias", "best_E_value", "best_score", "best_bias", "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", "description_of_target")
    read.table(filename, header = F, fill = T) %>%
        select(1:length(col_names)) %>%
        setNames(col_names)
}

pol <- read.table(pol_file) %>%
    select(Gene = 1, Domain = 2) %>%
    group_by(Gene) %>%
    summarize(has_exo = any(Domain == exo_domain), has_pol = any(Domain == pol_domain)) %>%
    mutate(Accession = sub("_[0-9]+$", "", Gene)) %>%
    group_by(Accession) %>%
    summarize(has_fused = any(has_pol & has_exo), has_pol = any(has_pol), has_exo = any(has_exo)) %>%
    mutate(Type = case_when(has_fused ~ "fused", has_pol & has_exo & !has_fused ~ "split", T ~ "unknown")) %>%
    with(setNames(Type, Accession))

psba <- read_tblout(psba_file) %>%
    mutate(Accession = sub("_[0-9]+$", "", target_name)) %>%
    distinct(Accession) %>%
    pull

strains <- read_excel(strain_file) %>%
    rename(Host = Strain, Host.taxonomy.prediction = Taxonomy)
clusters <- read_excel(cluster_file) %>%
    rename(Host.genus = Genus, Host.cluster = Cluster, Host.color = color)

col_names <- c("Accession", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
exo_nblas <- lapply(exo_nbla_files, read.table, header = F, row.names = NULL, col.names = col_names, sep = "\t") %>%
    bind_rows %>%
    filter(feature == "nblA") %>%
    select(Accession) %>%
    mutate(nblA = "exonuclease")

all_nblas <- read_jsonl(genome_nbla_file) %>%
    filter(full.score > nbla_min_score) %>%
    mutate(Accession = sub("_[0-9]+$", "", seq_name)) %>%
    filter(! Accession %in% exo_nblas$Accession) %>%
    mutate(nblA = "other") %>%
    bind_rows(exo_nblas)

ref_metadata <- read_excel(ref_file) %>%
    left_join(strains, by = "Host")
ingroup <- filter(ref_metadata, Group == "ingroup") %>%
    pull(Accession)
outgroup <- filter(ref_metadata, Group == "outgroup") %>%
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
add_mrca <- function(tree, colname, max_dist = 0) {
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
        mutate(mrca_col = ifelse(node == mrca & enough_tips, first(na.omit(my_column)), NA))
    tree[[paste0(colname, "_mrca")]] <- mrca$mrca_col
    mrca_data <- filter(mrca, !is.na(mrca_col))
    tip_data <- Descendants(phylo, mrca_data$node, type = "tips") %>%
        lapply(data.frame) %>%
        setNames(mrca_data$mrca_col) %>%
        bind_rows(.id = "mrca_col") %>%
        setNames(c("mrca_col", "node")) %>%
        mutate(label = phylo$tip.label[node]) %>%
        with(setNames(mrca_col, label))
    assigned <- mutate(tree, my_column = tip_data[label])
    tree[[paste0(colname, "_assigned")]] <- assigned$my_column
    return(tree)
}
annotate_tree <- function(tree, metadata, max_dist) {
    as_tibble(tree) %>%
        mutate(support = ifelse(node %in% parent, label, NA)) %>%
        separate_rows(support, sep = ";") %>%
        mutate(Accession = ifelse(node %in% parent, NA, label)) %>%
        separate(support, into = c("SH_aLRT", "UFBoot"), sep = "/", convert = T) %>%
        left_join(metadata, by = c(label = "Accession")) %>%
        arrange(-UFBoot) %>%
        distinct(node, .keep_all = T) %>%
        add_mrca(Clade, max_dist = max_dist) %>%
        to_treedata
}
to_treedata <- function(tree) {
    class(tree) <- c("tbl_tree", "tbl_df", "tbl", "data.frame")
    as.treedata(tree)
}   
keep_ingroup <- function(tree, ingroup) {
    tips_to_drop <- with(tree, tip.label[! tip.label %in% ingroup])
    tree <- drop.tip(tree, tips_to_drop, collapse.singles = F)
    singles <- Children(tree) %>%
        lapply(data.frame) %>%
        bind_rows(.id = "parent") %>%
        mutate(parent = as.integer(parent)) %>%
        setNames(c("parent", "node")) %>%
        group_by(parent) %>%
        filter(n() == 1, node > length(tree$tip.label))
    parents <- Ancestors(tree, type = "all") %>%
        lapply(data.frame) %>%
        bind_rows(.id = "node") %>%
        mutate(node = as.integer(node)) %>%
        setNames(c("node", "ancestor")) %>%
        filter(node %in% singles$node) %>%
        mutate(ancestor_is_not_single = ! ancestor %in% singles$parent) %>%
        filter(cumsum(ancestor_is_not_single) == 0) %>%
        mutate(parent_label = with(tree, c(tip.label, node.label))[ancestor]) %>%
        mutate(node_label   = with(tree, c(tip.label, node.label))[node]) %>%
        group_by(node) %>%
        summarize(label = paste(unique(c(node_label, parent_label)), collapse = ";"))
    tree$node.label[parents$node] <- parents$label
    collapse.singles(tree)
}
tree <- read.tree(tree_file) %>%
    root_tree(outgroup, ingroup)

imgvr_metadata <- read.csv(imgvr_file) %>%
    mutate(Accession = paste0(UVIG, "_1")) %>%
    mutate(Host.evidence = case_when(Topology == "Provirus" ~ "Prophage", Host.prediction.method == "" ~ NA_character_, T ~ Host.prediction.method))

gov2_metadata <- data.frame(Accession = tree$tip.label) %>%
    filter(grepl("^Station", Accession)) %>%
    mutate(Ecosystem.classification = "Environmental;Aquatic;Marine;")

env_metadata <- bind_rows(imgvr_metadata, gov2_metadata)

metadata <- bind_rows(ref_metadata, env_metadata) %>%
    complete(Accession = tree$tip.label) %>%
    distinct(Accession, .keep_all = T) %>%
    mutate(Ecosystem = sub(";[^;]*$", "", Ecosystem.classification)) %>%
    mutate(Ecosystem = ifelse(Ecosystem == ";;", NA, Ecosystem)) %>%
    mutate(Prophage = Topology == "Provirus") %>%
    get_taxon("phylum") %>%
    get_taxon("family") %>%
    get_taxon("genus") %>%
    left_join(clusters, by = "Host.genus") %>%
    left_join(all_nblas, by = "Accession") %>%
    distinct(Accession, Phage, Group, Clade, Host, Host.phylum, Host.genus, Host.cluster, Prophage, Host.evidence, nblA, Ecosystem) %>%
    mutate(psbA = Accession %in% psba) %>%
    mutate(polA = pol[Accession])

tree_dat_all <- tree %>%
    annotate_tree(metadata, max_mrca_dist)
tree_dat_ref <- keep_ingroup(tree, ingroup) %>%
    annotate_tree(metadata, max_mrca_dist)

write.table(as_tibble(tree_dat_all), file = "metadata.csv", sep = ",")

write.jtree(tree_dat_all, file = jtree_file)

plot_tree <- function(tree_dat, metadata, clusters) {
    evidence.shapes <- c(
        `Isolate taxonomy` = 15,
        `Prophage` = 15,
        `Related to isolate` = 0,
        `CRISPR spacer match` = 1,
        `k-mer match` = 1
    )
    cluster.colors <- with(clusters, setNames(Host.color, Host.cluster))
    ecosystems <- select(metadata, Accession, Ecosystem_tile = Ecosystem)
    ecosystem.colors <- c(
        `Environmental;Aquatic;Marine` = "blue",
        `Environmental;Aquatic;Freshwater` = "cyan"
    )
    polymerase <- select(metadata, Accession, polA_tile = polA) %>%
        filter(!is.na(polA_tile))
    ggtree(tree_dat, size = 0.1, layout = "rectangular", open.angle = 45) +
        geom_nodepoint(aes(x = branch, subset = !is.na(UFBoot) & UFBoot >= 95), size = 0.2, color = "#4d4d4dff") +
        geom_tiplab(aes(subset = !is.na(Phage), label = Phage, color = Host.cluster), size = 2, linesize = 0.1, fontface = "bold", align = T) +
        scale_color_manual(values = cluster.colors) + new_scale_color() +
        geom_tippoint(aes(subset = !is.na(nblA), color = nblA), size = 0.5) + new_scale_color() +
        geom_fruit(aes(x = 0, y = Accession, fill = Ecosystem_tile), ecosystems, geom_tile, axis.params = c(axis = "x", text = "Ecosystem"), offset = 0.5, width = 0.1) +
        scale_fill_manual(values = ecosystem.colors) + new_scale_fill() +
        geom_tippoint(aes(subset = !is.na(psbA) & psbA, x = x + 0.02), color = "#009845", shape = 15, size = 0.5) +
        geom_fruit(aes(x = 0, y = Accession, fill = polA_tile), polymerase, geom_tile, axis.params = c(axis = "x", text = "Polymerase"), offset = 0.5, width = 0.1) +
        geom_cladelab(mapping = aes(subset = !is.na(Clade_mrca), node = node, label = Clade_mrca), align = T, offset = 0.3) +
        geom_treescale(width = 0.4)
}

p_all <- plot_tree(tree_dat_all, metadata, clusters)
p_ref <- plot_tree(tree_dat_ref, metadata, clusters)

ggsave(plot_all_file, p_all, width = w_all, height = h_all)
ggsave(plot_ref_file, p_ref, width = w_ref, height = h_ref)
