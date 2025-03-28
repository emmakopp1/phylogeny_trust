#!/bin/sh
# Fichier sim-17.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-17.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-17/ctmc-strict-bd-17.xml