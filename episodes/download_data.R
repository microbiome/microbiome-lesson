# Helper used when the lesson is built.
#
# Each episode sources this file in a hidden chunk so that it can be built on
# its own. It recreates the objects that earlier episodes save to disk.
# Learners do not need this file: they create the same objects themselves by
# following the episodes in order.

suppressPackageStartupMessages(library(mia))

# Episode 04: import the data into a TreeSE
if (!file.exists(file.path("data", "tse.rds"))) {
    counts  <- read.csv(file.path("data", "taxonomy_abundance.csv"), row.names = 1)
    taxa    <- read.csv(file.path("data", "taxonomy_table.csv"), row.names = 1)
    samples <- read.csv(file.path("data", "sample_metadata.csv"), row.names = 1)
    tree    <- ape::read.tree(file.path("data", "phylogeny.tree"))
    tse <- TreeSummarizedExperiment(
        assays  = list(counts = as.matrix(counts)),
        rowData = taxa,
        colData = samples,
        rowTree = tree,
        mainExpName = "main"
    )
    saveRDS(tse, file.path("data", "tse.rds"))
}

# Episode 05: add library size and transformed tables
if (!file.exists(file.path("data", "tse_processed.rds"))) {
    tse <- readRDS(file.path("data", "tse.rds"))
    tse$library_size <- colSums(assay(tse, "counts"))
    tse <- transformAssay(tse, assay.type = "counts", method = "relabundance")
    tse <- transformAssay(tse, assay.type = "counts", method = "rclr")
    saveRDS(tse, file.path("data", "tse_processed.rds"))
}
