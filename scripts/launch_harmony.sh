#!/usr/bin/env bash
### PBS HEADER
#PBS -N harmony
#PBS -o logs/harmony.log
#PBS -j oe

#PBS -q batch
#PBS -l walltime=00:05:00
#PBS -l nodes=1:ppn=1
### PBS HEADER

cd $PBS_O_WORKDIR

module load R/3.5.1

echo $INPUT_DIR
echo $OUTPUT_DIR
Rscript scripts/run_harmony.R $INPUT_DIR $OUTPUT_DIR
