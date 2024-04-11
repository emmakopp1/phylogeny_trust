#!/bin/sh
# Fichier ex.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=ex.out

~/beast/bin/beast -overwrite /home/users/kopp/work/chr_paper/reconstruction/one-tree/beast-data-sim-8/ctmc-strict-bd-8.xml