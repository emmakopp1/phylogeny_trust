#!/bin/bash
#SBATCH --job-name "ancr.sbatch"
#SBATCH --cpus-per-task=50
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --partition=debug
#SBATCH --output=/home/users/kopp/work/ancestral_reconstruction/test.out
#SBATCH --error=/home/users/kopp/work/ancestral_reconstruction/test.out

Rscript /home/users/kopp/work/ancestral_reconstruction/ancr_limit_age_parralel.R