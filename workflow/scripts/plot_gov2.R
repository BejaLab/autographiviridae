library <- function(...) suppressPackageStartupMessages(base::library(...))

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(grid)
library(scatterpie)
library(tools)
library(RColorBrewer)

with(snakemake@input, {
    txt_files <<- txt
    stat_files <<- stat
    metadata_file <<- metadata
})
with(snakemake@output, {
    srf_file <<- srf
})
with(snakemake@params, {
    sras <<- sras
})

stats <- lapply(stat_files, read.table, header = T) %>%
    setNames(sras) %>%
    bind_rows(.id = "sra") %>%
    mutate(num_seqs = as.numeric(gsub(",", "", num_seqs)))

metadata <- read_excel(metadata_file, skip = 1) %>%
    extract(`Environmental feature  (http://environmentontology.org/)`, into = "Layer", regex = "\\[([A-Z]{3})\\]") %>%
    rename(sra = `Run accession number(s) registered at the European Nucleotides Archive (www.ebi.ac.uk/ena)`) %>%
    separate_rows(sra) %>%
    mutate(depth = as.numeric(`Depth (m)`)) %>%
    mutate(x = as.numeric(`Longitude (decimal degree E)`), y = as.numeric(`Latitude (decimal degree N)`))

data <- lapply(txt_files, read.table, header = T) %>%
    lapply(select, c(Geneid, Chr, Start, End, Length, Count = 7)) %>%
    setNames(sras) %>%
    bind_rows(.id = "sra") %>%
    separate_rows(Start, End, Chr, sep = ";", convert = T) %>%
    group_by(Geneid, sra, Count) %>%
    summarize(Median_len = median(End - Start + 1), .groups = "drop") %>%
    left_join(stats, by = "sra") %>%
    left_join(metadata, by = "sra") %>%
    group_by(Station, Layer, Geneid, x, y, depth) %>%
    summarize(FPKM = 10^9 * sum(Count) / mean(Median_len) / sum(num_seqs), .groups = "drop") %>%
    separate(Geneid, into = c("Gene", "Clade"), sep = "_") %>%
    spread(Gene, FPKM, fill = 0) %>%
    mutate(Total = pmax(exonuclease, nblA), no_nblA = Total - nblA) %>%
    gather(Category, FPKM, nblA, no_nblA) %>%
    mutate(Clade_Category = paste(Clade, Category, sep = "@")) %>%
    select(Station, Layer, Clade_Category, x, y, depth, FPKM) %>%
    spread(Clade_Category, FPKM) %>%
    mutate(Total = 1.5*log10(rowSums(across(contains("@"))))) %>%
    mutate(Total = ifelse(Total == -Inf, 0, Total))
data_layers <- split(data, f = data$Layer)

clade_categories <- names(data)[grepl("@", names(data))] %>%
    setNames(brewer.pal(length(.), "Paired"),.)

get_repel_coords <- function(.data, map_g, width, height) {
    grid.newpage()
    pushViewport(viewport(width = width, height = height))
    g <- map_g +
        geom_point(aes(x, y), data = .data) +
        geom_text_repel(aes(x, y, size = Total), size = 3, label = ".", data = .data, box.padding = 5, max.overlaps = Inf)
    panel_params <- ggplot_build(g)$layout$panel_params[[1]]
    xrg <- panel_params$x.range
    yrg <- panel_params$y.range

    textrepeltree <- ggplotGrob(g) %>%
        grid.force(draw = F) %>%
        getGrob("textrepeltree", grep = T)
    children <- childNames(textrepeltree) %>%
        grep("textrepelgrob", ., value = T)

    get_xy <- function(n) {
        grob <- getGrob(textrepeltree, n)
        data.frame(
            x.repel = xrg[1] + diff(xrg) * convertX(grob$x, "native", valueOnly = T),
            y.repel = yrg[1] + diff(yrg) * convertY(grob$y, "native", valueOnly = T)
        )
    }
    lapply(children, get_xy) %>%
        bind_rows %>%
        cbind(.data) %>%
        mutate(theta = atan2(y - y.repel, x - x.repel), x.segm = x.repel + Total * cos(theta), y.segm = y.repel + Total * sin(theta))
}
my_ggplot_world <- function() {
    world <- map_data("world")
    ggplot(world) +
        geom_map(data = world, map = world, aes(map_id = region), color = "black", fill = "lightgray", linewidth = 0.1) +
        coord_equal() +
        xlim(c(-185, 185))
}

g <- my_ggplot_world()

width <- 16
height <- 8

data_repel <- get_repel_coords(data_layers$SRF, g, width, height)
p <- g +
    geom_segment(aes(x = x.segm, y = y.segm, xend = x, yend = y), data = data_repel, arrow = arrow(length = unit(0.01, "npc"))) +
    geom_scatterpie(aes(x = x.repel, y = y.repel, r = Total), data = data_repel, color = NA, cols = names(clade_categories)) +
    # geom_circle(aes(x0 = x.repel, y0 = y.repel, r = 0.01), data = data_repel) +
    # geom_scatterpie_legend(.data$Total, x = -140, y = -70, labeller = labeller) +
    scale_fill_manual(values = clade_categories) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.title.y = element_blank())
ggsave(srf_file, p, width = width, height = height)
