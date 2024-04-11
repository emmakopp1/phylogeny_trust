#!/bin/sh
# Fichier sino-tibet.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=ex.out

~/beast/bin/beast -overwrite /home/users/kopp/indo-euro/chr_iecor/IECoR2-chr.xml