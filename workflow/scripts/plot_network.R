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
    metadata_file <<- metadata
    network_file <<- network   # analysis/nbla/nbla_all_cdhit.nex
    clstr_file <<- clstr # analysis/nbla/nblA_all_cdhit.faa.clstr
})
with(snakemake@output, {
    plot_file <<- plot
    data_file <<- data
})

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
metadata <- read.csv(metadata_file, na.strings = "") %>%
    rename(Group = nbla_category, Category = taxon)
network <- read.nexus.networx(network_file)

node_data <- read_clstr(clstr_file) %>%
    mutate(Accession = sub("_[0-9]+$", "", label)) %>%
    mutate(Prefix = sub("^([^_]+_[^_]+).+", "\\1", Accession)) %>%
    left_join(metadata, by = "label") %>%
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
write.csv(node_data, data_file)
p <- ggsplitnet(network, size = 0.1) %>% add_node_data(node_data) +
    geom_point2(aes(subset = !is.na(Group), x = x, y = y, color = Group, shape = Category)) +
    geom_tiplab2(aes(label = alias, x = x_ext(xend, yend, x, y, 0.03), y = y_ext(xend, yend, x, y, 0.03)), color = "gray30") +
    new_scale_color() +
    geom_treescale(width = 0.1) +
    coord_fixed()
ggsave(plot_file, p, width = 10, height = 8)
