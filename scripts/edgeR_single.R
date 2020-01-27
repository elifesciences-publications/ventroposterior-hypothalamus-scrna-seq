#!/usr/bin/env Rscript
#.libPaths(.libPaths()[2])
setwd("/projects/robson-lab/research/hypothalamus/")

library("edgeR")
library("dplyr")

args = commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
group1 <- args[2]
group2 <- args[3]

raw.counts <- t(
  read.csv(
    paste(input_dir, "counts.csv", sep="/"),
    row.names = NULL,
    header = FALSE
  )
)

obs <- read.csv(paste(input_dir, "obs.csv", sep="/"), row.names=1)
cluster.inds <- obs[, "cluster_revised"] %in% c(group1, group2)
cluster.cells <- rownames(obs[cluster.inds, ])

var <- read.csv(paste(input_dir, "var.csv", sep="/"), row.names = 1)

colnames(raw.counts) <- rownames(obs)
raw.counts <- as.data.frame(raw.counts,
                            row.names = rownames(var))
rownames(raw.counts) <- rownames(var)
obs <- obs[cluster.inds, ]
raw.counts <- raw.counts[, cluster.cells]

datestamp <- format(Sys.time(), "%Y-%m-%dT%H-%M")
fname <- paste0("dge_neuronal_", group1, "-v-", group2, "_", datestamp, ".csv")
outfile <- paste(input_dir, fname, sep="/")

group.ids <- obs[, "cluster_revised"]

group <- factor(group.ids, levels = c(group1, group2))
print(length(group.ids))
print(group)
design <- model.matrix(~group)
print(design)

cluster.counts <- raw.counts
print(dim(cluster.counts))
keep <- rowSums(cluster.counts > 1) >= 2
cluster.counts <- cluster.counts[keep, ]
print(dim(cluster.counts))

y <- DGEList(counts = cluster.counts, group = group)
print("Made DGEList")

y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y)
y <- estimateGLMTagwiseDisp(y)
print("Got dispersions")

fit <- glmFit(y, design)
lrt <- glmLRT(fit)#, coef=2)
print("Completed fit")

dge <- topTags(lrt, n = nrow(y$counts))
dge <- dge$table
print("Writing DGE list top tags")

write.csv(dge, file=outfile)
