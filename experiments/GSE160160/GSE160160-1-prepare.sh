#!/bin/bash

#SBATCH --job-name=GSE160160-prepare
#SBATCH --output=GSE160160-1-prepare.out
#SBATCH --error=GSE160160-1-prepare.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500G
#SBATCH --partition=owners

ml R/4.2.0
ml biology physics geos/3.6.2 hdf5/1.12.2
Rscript GSE160160-1-prepare.R
