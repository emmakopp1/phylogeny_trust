#!/bin/sh
# Fichier bantu-subsample.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=%x%J.out

~/beast/bin/beast -overwrite /home/users/kopp/work/bantu/bantu-ctmc-strict-bd-subsample/bantu-ctmc-strict-bd-subsample.xml