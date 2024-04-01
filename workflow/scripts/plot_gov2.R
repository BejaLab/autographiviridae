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
library(jsonlite)
library(ggforce)
library(rgdal)

if (interactive()) {
    txt_files <- txt_files <- Sys.glob("analysis/recruitment/GOV2/bwa/ERR*.txt")
    sras <- sub(".txt", "", basename(txt_files))
    stat_files <- Sys.glob("fastq/GOV2/*.stat")
    metadata_file <- "analysis/recruitment/GOV2/metadata.xlsx"
    srf_file <- "SRf.pdf"
    dcm_file <- "DCM.pdf"
    pie_file <- "pies.pdf"
    jtree_file <- "output/phylophlan_markers9.jtree"
    clade_colors <- list(
      A = c( '#33a02cff', '#b2df8aff' ),
      B = c( '#0000ffff', '#96c3dcff' ),
      C = c( '#d90017ff', '#f88587ff' ),
      D = c( '#fc6908ff', '#fbb25cff' ),
      E = c( '#d45e9eff', '#eca4d1ff' ),
      F = c( '#562888ff', '#bea0ccff' )
    )
} else {
    with(snakemake@input, {
        txt_files <<- txt
        stat_files <<- stat
        metadata_file <<- metadata
        jtree_file <<- jtree
        shp_file <<- shp
    })
    with(snakemake@output, {
        srf_file <<- srf
        dcm_file <<- dcm
        pie_file <<- pies
        pie_data_file <<- pie_data
    })
    with(snakemake@params, {
        sras <<- sras
        clade_colors <<- colors
    })
}

forward_scale <- function(x, offset = 10) {
    ifelse(x == 0, 0, log(x + offset))
}
reverse_scale <- function(x, offset = 10, digits = 5) {
    ifelse(x == 0, 0, round(exp(x) - offset, digits))
}

genome_data <- fromJSON(jtree_file)

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
    rename(nblA_plus = nblA) %>%
    mutate(Total = pmax(exonuclease, nblA_plus), nblA_minus = Total - nblA_plus) %>%
    gather(Category, FPKM, nblA_plus, nblA_minus) %>%
    mutate(Clade_Category = paste(Clade, Category, sep = "@")) %>%
    select(Station, Layer, Clade_Category, x, y, depth, FPKM) %>%
    spread(Clade_Category, FPKM) %>%
    mutate(Total = rowSums(across(contains("@")))) %>%
    mutate(nblA_plus = rowSums(across(contains("@nblA_plus")))) %>%
    mutate(Size = forward_scale(Total)) %>%
    mutate(nblA_pct = nblA_plus / Total)

data_layers <- split(data, f = data$Layer)

clade_categories <- lapply(clade_colors, t) %>%
    lapply(data.frame) %>%
    bind_rows(.id = "clade") %>%
    setNames(c("clade", "nblA_plus", "nblA_minus")) %>%
    gather(type, color, -clade) %>%
    mutate(category = paste(clade, type, sep = "@")) %>%
    with(setNames(color, category))
get_repel_coords <- function(.data, map_g, width, height) {
    grid.newpage()
    pushViewport(viewport(width = width, height = height))
    g <- map_g +
        geom_text_repel(aes(x, y, size = Size*22, point.size = Size), label = "O", data = .data, box.padding = 0, max.overlaps = Inf) +
        scale_size_identity()
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
        mutate(theta = atan2(y - y.repel, x - x.repel), x.segm = x.repel + Size * cos(theta), y.segm = y.repel + Size * sin(theta))
}
my_ggplot_world <- function(shape_file) {
    map_base <- file_path_sans_ext(shape_file)
    world.land <- readOGR(dirname(shape_file), basename(file_path_sans_ext(shape_file)))
    continents.simple <- fortify(world.land)
    all.continents <- data_frame(id = unique(as.character(continents.simple$id)))
    ggplot(all.continents) +
        geom_map(data = continents.simple, map = continents.simple, aes(map_id = id), color = "lightgray", fill = "lightgray") +
        coord_quickmap(xlim = c(-185, 185), ylim = c(-85, 90), expand = F)
}
g <- my_ggplot_world(shp_file)

size_breaks <- c(10, 100, 1000, 10000)
file_names <- list(
    SRF = srf_file,
    DCM = dcm_file
)

genome_pie <- genome_data$data %>%
    filter(!is.na(Accession)) %>%
    replace_na(list(Clade_assigned = "unknown")) %>%
    mutate(Category = recode(nblA, exonuclease = "nblA_plus", .missing = "nblA_minus")) %>%
    mutate(Clade_Category = paste(Clade_assigned, Category, sep = "@")) %>%
    group_by(Clade_Category) %>%
    summarize(n = n()) %>%
    mutate(pct = 100 * n / sum(n), Layer = "Whole genomes")
fpkm_pie <- data_layers %>%
    lapply(select, all_of(names(clade_categories))) %>%
    lapply(colSums) %>%
    bind_rows(.id = "Layer") %>%
    gather(Clade_Category, FPKM, -Layer) %>%
    group_by(Layer) %>%
    mutate(pct = 100 * FPKM / sum(FPKM)) %>%
    filter(Layer %in% names(file_names))
pie_data <- bind_rows(genome_pie, fpkm_pie)

# Get the positions
pie_labels <- group_by(pie_data, Layer) %>%
    filter(pct > 0) %>%
    mutate(csum = rev(cumsum(rev(pct))), pos = pct/2 + lead(csum, 1), pos = ifelse(is.na(pos), pct/2, pos))

pie <- ggplot(pie_data, aes(x = "", y = pct, fill = Clade_Category)) +
    geom_bar(stat = "identity", width = 1) +
    geom_text(aes(label = paste0(round(pct, 2), "%")), size = 3, position = position_stack(vjust = 0.5)) +
    scale_fill_manual(values = clade_categories) +
    coord_polar("y", start = 0) +
    facet_wrap(~ Layer, ncol = 2) +
    geom_label_repel(data = pie_labels, aes(label = Clade_Category, y = pos, x = after_stat(1.33)),
        point.padding = NA,
        max.overlaps = Inf,
        nudge_x = 0.5,
        force = 1,
        force_pull = 0,
        segment.alpha = 0.5
    ) +
    theme_void() +
    theme(legend.position = "None")
ggsave(pie_file, pie)

write.csv(pie_data, file = pie_data_file)

invisible(
    lapply(names(file_names), function(layer) {
        width <- 16
        height <- 8
        data_repel <- get_repel_coords(data_layers[[layer]], g, width, height)
        p <- g +
            geom_segment(aes(x = x.segm, y = y.segm, xend = x, yend = y), data = data_repel, arrow = arrow(length = unit(0.01, "npc"))) +
            geom_scatterpie(aes(x = x.repel, y = y.repel, r = Size), data = data_repel, color = NA, cols = names(clade_categories)) +
            geom_circle(aes(x0 = x.repel, y0 = y.repel, r = Size), data = data_repel) +
            geom_scatterpie_legend(data$Size, x = -140, y = -70, breaks = forward_scale(size_breaks), labeller = reverse_scale) +
            scale_fill_manual(values = clade_categories) +
            theme_bw() +
            theme(axis.title.x = element_blank(), axis.title.y = element_blank())
        ggsave(file_names[[layer]], p, width = width, height = height)
    })
)

#data_depth <- filter(data, !is.na(depth), depth <= 200, Total > 0)
#p <- ggplot(data_depth) +
#    geom_point(aes(x = nblA_pct * 100, y = depth, size = Total), data = data_depth) +
#    #geom_scatterpie(aes(x = nblA_pct * 100, y = depth, r = Size), data = data, color = NA, cols = names(clade_categories)) +
#    #geom_circle(aes(x0 = nblA_pct * 100, y0 = depth, r = Size), data = data) +
#    #geom_scatterpie_legend(data$Size, x = 0.5, y = 0, breaks = forward_scale(size_breaks), labeller = reverse_scale) +
#    scale_y_continuous(trans = trans_reverser('log10')) +
#    scale_fill_manual(values = clade_categories) +
#    theme_bw()
#    #coord_fixed()
#ggsave("tmp.pdf", p, width = 5, height = 10)
