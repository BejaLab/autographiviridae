library <- function(...) suppressPackageStartupMessages(base::library(...))
library(DiffLogo)
library(treeio)
library(dplyr)
library(tidyr)

with(snakemake@input, {
    tree_file <<- tree
    aln_file <<- aln
})
with(snakemake@params, {
    group1_clades <<- group1
    group2_clades <<- group2
    outliers <<- outliers
    outgroup_prefix <<- outgroup
})
output_file <- unlist(snakemake@output)

zappo <- c(
    ILVAM = "#ffafaf",
    FWY = "#ffc800",
    KRH = "#6464ff",
    DE = "#ff0000",
    STNQ = "#00ff00",
    PG = "#ff00ff",
    C = "#ffff00"
)
alphabet <- data.frame(chars = names(zappo), cols = unname(zappo)) %>%
    separate_rows(chars, sep = "") %>%
    filter(chars != "") %>%
    arrange(chars) %>%
    as.list
alphabet$supportReverseComplement <- F
alphabet$size <- length(alphabet$chars)
class(alphabet) <- "Alphabet"

aln <- read.fasta(aln_file, type = "AA") %>%
    as.character %>%
    lapply(paste, collapse = "") %>%
    unlist %>%
    toupper
names(aln) <- sub(" .*", "", names(aln))
aln <- aln[!grepl("_consensus", names(aln))]

tree <- as_tibble(read.jtree(tree_file)) 

to_pwm <- function(labels, aln, alphabet) {
    labels <- labels[labels %in% names(aln)]
    getPwmFromAlignment(aln[labels], alphabet = alphabet) %>%
        as.data.frame
}
select_labels <- function(tree, clades, outliers) {
    filter(tree, Clade_assigned %in% clades, ! label %in% outliers) %>%
        pull(label)
}
get_max <- function(pwm) {
    rownames_to_column(pwm, "aa") %>%
        gather(pos, p, -aa) %>%
        arrange(-p) %>%
        group_by(pos) %>%
        summarize(p = first(p)) %>%
        mutate(pos = as.integer(pos)) %>%
        arrange(pos) %>%
        pull(p)
}

group1_labels <- select_labels(tree, group1_clades, outliers)
group1 <- to_pwm(group1_labels, aln, alphabet)
group2_labels <- select_labels(tree, group2_clades, outliers)
group2 <- to_pwm(group2_labels, aln, alphabet)
outgroup_labels <- names(aln) %>%
    `[`(grepl(outgroup_prefix, .))
outgroup <- to_pwm(outgroup_labels, aln, alphabet)

diffObj <- prepareDiffLogoTable(list(group1 = group1, group2 = group2), alphabet = alphabet, configuration = diffLogoTableConfiguration(alphabet = alphabet, enableClustering = F))
diffObj$diffLogoObjMatrix <- enrichDiffLogoTableWithPvalues(diffObj$diffLogoObjMatrix, c(group1 = length(group1_labels), group2 = length(group2_labels)))

data <- data.frame(group1 = get_max(group1), group2 = get_max(group2)) %>%
    mutate(pvals12 = diffObj$diffLogoObjMatrix$group1$group2$pvals) %>%
    mutate(pvals21 = diffObj$diffLogoObjMatrix$group2$group1$pvals) %>%
    mutate(js = diffObj$diffLogoObjMatrix$group2$group1$heights) %>%
    mutate(pos = 1:n())

pwms <- list(
    group1 = group1,
    group2 = group2,
    outgroup = outgroup
)
sizes <- c(
    group1 = length(group1_labels),
    group2 = length(group2_labels),
    outgroup = length(outgroup_labels)
)

ss <- data.frame(pred = aln["ss_pred"], conf = aln["ss_conf"]) %>%
    separate_rows(pred, conf, sep = "") %>%
    filter(pred != "") %>%
    mutate(conf = as.integer(conf) / 9) %>%
    mutate(pos = 1:n()) %>%
    spread(pred, conf, fill = 0) %>%
    select(-pos) %>%
    mutate(Y = 1 - rowSums(.))
missing <- setdiff(alphabet$chars, names(ss))
ss[missing] <- 0
ss <- ss[ ,order(names(ss))]  %>%
    t %>% as.data.frame %>%
    `names<-`(1:ncol(.))

pdf(output_file)
seqLogo(ss, alphabet = alphabet, stackHeight = function(x) list(height = x * 4.5, ylab = ""))
diffLogoTable(PWMs = pwms, alphabet = alphabet, configuration = diffLogoTableConfiguration(alphabet = alphabet, enableClustering = F))
drawDiffLogoTable(diffObj)
dev.off()
