#!/bin/sh
# Fichier trans.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=trans.out

~/beast/bin/beast -overwrite /home/users/kopp/work/trans-ctmc-strict-fbd/ctmc-strict-fbd.xml