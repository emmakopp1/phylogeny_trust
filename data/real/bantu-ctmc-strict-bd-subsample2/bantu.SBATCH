#!/bin/sh
# Fichier bantu-subsample2.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=%x%J.out

~/beast/bin/beast -overwrite /home/users/kopp/work/bantu/bantu-ctmc-strict-bd-subsample2/bantu-ctmc-strict-bd-subsample2.xml