library <- function(...) suppressPackageStartupMessages(base::library(...))
library(DiffLogo)
library(dplyr)
library(tidyr)
library(Biostrings)

with(snakemake@input, {
    metadata_file <<- metadata
    aln_file <<- aln
})
with(snakemake@params, {
    min_seqs <<- min_seqs
    js_threshold <<- 0.9
})
output_file <- unlist(snakemake@output)

groups <- read.csv(metadata_file) %>%
    select(label, category = nbla_category2) %>%
    group_by(category) %>%
    filter(n() >= min_seqs) %>%
    split(f = .$category) %>%
    lapply(pull, "label")

zappo <- c(
    ILVAM = "#ffafaf",
    FWY = "#ffc800",
    KRH = "#6464ff",
    DE = "#ff0000",
    STNQ = "#00ff00",
    PG = "#ff00ff",
    C = "#ffff00",
    X = "black"
)
alphabet <- data.frame(chars = names(zappo), cols = unname(zappo)) %>%
    separate_rows(chars, sep = "") %>%
    filter(chars != "") %>%
    arrange(chars) %>%
    as.list
alphabet$supportReverseComplement <- F
alphabet$size <- length(alphabet$chars)
class(alphabet) <- "Alphabet"

aln <- readBStringSet(aln_file) %>%
    as.character %>%
    lapply(function(x) gsub("[a-z]", "", x)) %>%
    lapply(function(x) gsub("-", "X", x))
names(aln) <- sub(" .*", "", names(aln))

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
labels <- lapply(names(groups), function(g) {
    groups[[g]][groups[[g]] %in% names(aln)]
}) %>% setNames(names(groups))
pwms <- lapply(labels, function(group_labels) {
    pwm <- getPwmFromAlignment(aln[group_labels], alphabet = alphabet)
    as.data.frame(pwm)
})
sizes <- unlist(lapply(labels, length))

diffObj <- prepareDiffLogoTable(pwms, alphabet = alphabet, configuration = diffLogoTableConfiguration(alphabet = alphabet, enableClustering = T))
diffObj$diffLogoObjMatrix <- enrichDiffLogoTableWithPvalues(diffObj$diffLogoObjMatrix, sizes)
# trick drawDiffLogoTable into showing asterisks for positions with JS divergence above threshold
diffObj$diffLogoObjMatrix <- lapply(diffObj$diffLogoObjMatrix, function(x) lapply(x, function(y) {
    y$pvals <- ifelse(y$heights > js_threshold & y$pvals < 0.01, 0, 1)
    return(y)
}))

ss <- data.frame(pred = aln$ss_pred, conf = aln$ss_conf) %>%
    separate_rows(pred, conf, sep = "") %>%
    filter(pred != "") %>%
    mutate(conf = as.integer(conf) / 9) %>%
    mutate(pos = 1:n()) %>%
    spread(pred, conf, fill = 0) %>%
    select(-pos) %>%
    mutate(Y = 1 - rowSums(.)) # dummy character to guarantee sum == 1
missing <- setdiff(alphabet$chars, names(ss))
ss[missing] <- 0
ss <- ss[ ,order(names(ss))]  %>%
    t %>% as.data.frame %>%
    `names<-`(1:ncol(.))

pdf(output_file)
seqLogo(ss, alphabet = alphabet, stackHeight = function(x) list(height = x * 4.5, ylab = ""))
# diffLogoTable(PWMs = pwms, alphabet = alphabet, configuration = diffLogoTableConfiguration(alphabet = alphabet, enableClustering = T))
drawDiffLogoTable(diffObj)
dev.off()
