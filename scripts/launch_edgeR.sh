#!/usr/bin/env bash
#### PBS preamble

#PBS -N edgeR
#PBS -o /projects/robson-lab/research/hypothalamus/logs/log.edgeR.${PBS_JOBID%%.*}.out

#PBS -j oe
#PBS -m a

#PBS -q batch
#PBS -l vmem=16GB
#PBS -l nodes=1:ppn=1
#PBS -l walltime=01:00:00

#### PBS preamble

cd /projects/robson-lab/research/hypothalamus/

module load R/3.5.1

Rscript scripts/edgeR_single.R final_analysis/neuronal_13-17_dge 1 2
Rscript scripts/edgeR_single.R final_analysis/neuronal_13-17_dge 1 3
Rscript scripts/edgeR_single.R final_analysis/neuronal_13-17_dge 1 4
Rscript scripts/edgeR_single.R final_analysis/neuronal_13-17_dge 2 3
Rscript scripts/edgeR_single.R final_analysis/neuronal_13-17_dge 2 4
Rscript scripts/edgeR_single.R final_analysis/neuronal_13-17_dge 3 4
