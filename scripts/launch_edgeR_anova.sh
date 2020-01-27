#!/usr/bin/env bash
#### PBS preamble

#PBS -N edgeR
#PBS -o /projects/robson-lab/research/hypothalamus/logs/log.edgeR-anova.${PBS_JOBID%%.*}.out

#PBS -j oe
#PBS -m a

#PBS -q batch
#PBS -l vmem=24GB
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00

#### PBS preamble

cd /projects/robson-lab/research/hypothalamus

module load R/3.5.1

#Rscript scripts/edgeR_anova.R final_analysis/ion_channel "ion"
#Rscript scripts/edgeR_anova.R final_analysis/transcription_factor "tfs"
#Rscript scripts/edgeR_anova.R final_analysis/dnabinding_transcription_factor "dbtfs"
Rscript scripts/edgeR_anova.R final_analysis/neuronal_13-17_dge "13-vs-17"
