#!/usr/bin/env Rscript
library(harmony)

setwd("/projects/robson-lab/research/hypothalamus/final_analysis")

args = commandArgs(trailingOnly=TRUE)
input_dir <- args[1]
output_dir <- args[2]

pc.embeddings <- read.csv(
    paste(input_dir, "pc_embeddings.csv", sep="/"),
    header=FALSE,
)
metadata <- read.csv(
    paste(input_dir, "metadata.csv", sep="/"),
    header=TRUE,
    stringsAsFactors=FALSE
)
covariates <- read.csv(
    paste(input_dir, "covariates.csv", sep="/"),
    header=FALSE,
    stringsAsFactors=FALSE
)
thetas <- read.csv(
    paste(input_dir, "thetas.csv", sep="/"),
    header=FALSE,
    stringsAsFactors=FALSE
)
covariates <- as.vector(t(covariates))
thetas <- as.vector(t(thetas))
print(covariates)
print(thetas)

harmony_embeddings <- HarmonyMatrix(
    pc.embeddings,
    metadata,
    vars_use=covariates,
    theta=thetas,
    do_pca=FALSE
)

write.csv(
    harmony_embeddings,
    file=paste(output_dir, "harmony_pc_embeddings.csv", sep="/")
)
