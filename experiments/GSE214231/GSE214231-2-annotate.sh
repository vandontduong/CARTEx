#!/bin/bash

#SBATCH --job-name=GSE214231-annotate
#SBATCH --output=GSE214231-2-annotate.out
#SBATCH --error=GSE214231-2-annotate.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500G
#SBATCH --partition=owners

ml R/4.2.0
ml biology physics geos/3.6.2 hdf5/1.12.2
Rscript GSE214231-2-annotate.R
