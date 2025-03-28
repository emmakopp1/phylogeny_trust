#!/bin/sh
# Fichier sim-16.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-16.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-16/ctmc-strict-bd-16.xml