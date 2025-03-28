#!/bin/sh
# Fichier sim-8.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-8.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-8/ctmc-strict-bd-8.xml