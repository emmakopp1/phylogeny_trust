#!/bin/sh
# Fichier sim-7.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-7.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-7/ctmc-strict-bd-7.xml