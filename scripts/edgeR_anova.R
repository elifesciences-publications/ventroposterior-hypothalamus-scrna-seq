#!/usr/bin/env Rscript
#.libPaths(.libPaths()[2])
setwd("/projects/robson-lab/research/hypothalamus")

library("edgeR")
library("dplyr")

args = commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
out_name <- args[2]

counts_file <- paste(input_dir, "counts.csv", sep="/")
obs_file <- paste(input_dir, "obs.csv", sep="/")
var_file <- paste(input_dir, "var.csv", sep="/")

datestamp <- format(Sys.time(), "%Y-%m-%dT%H-%M")
fname <- paste0("dge_anova_", out_name, "_", datestamp, ".csv")
outfile <- paste(input_dir, fname, sep="/")

# Load counts
raw.counts <- t(read.csv(counts_file, row.names = NULL, header = FALSE))
print(dim(raw.counts))

# Load obs and var
obs <- read.csv(obs_file, row.names = 1)
var <- read.csv(var_file, row.names = 1)

# Create counts DF
colnames(raw.counts) <- rownames(obs)
raw.counts <- as.data.frame(raw.counts,
                            row.names = rownames(var))
rownames(raw.counts) <- rownames(var)

# Define groups
group.ids <- obs[, "cluster_revised"]
groups <- unique(group.ids)
group <- factor(group.ids, levels = groups)
design <- model.matrix(~group)

# Do some QC
counts <- raw.counts
keep <- rowSums(counts > 1) >= 2
counts <- counts[keep, ]

y <- DGEList(counts = counts, group = group)
print("Made DGEList")

y <- calcNormFactors(y)
y <- estimateGLMCommonDisp(y)
y <- estimateGLMTagwiseDisp(y)
print("Got dispersions")

#fit <- glmFit(y, design)
#lrt <- glmLRT(fit)#, coef=2)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2:length(groups))
print("Completed fit")

#dge <- topTags(lrt, n = nrow(y$counts))
dge <- topTags(qlf, n = nrow(y$counts))
dge <- dge$table
print("Writing DGE list top tags")

write.csv(dge, file=outfile)
