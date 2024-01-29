#!/bin/bash

#SBATCH --job-name=GSE136184-l-prepare
#SBATCH --output=GSE136184-l-1-prepare.out
#SBATCH --error=GSE136184-l-1-prepare.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500G
#SBATCH --partition=owners

ml R/4.2.0
ml biology physics geos/3.6.2 hdf5/1.12.2
Rscript GSE136184-l-1-prepare.R
