#!/bin/bash

#SBATCH --job-name=Zenodo3993994-csf-annotate
#SBATCH --output=Zenodo3993994-csf-2-annotate.out
#SBATCH --error=Zenodo3993994-csf-2-annotate.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=500G
#SBATCH --partition=owners

ml R/4.2.0
ml biology physics geos/3.6.2 hdf5/1.12.2
Rscript Zenodo3993994-csf-2-annotate.R
