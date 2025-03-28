#!/bin/sh
# Fichier sim-1.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-1.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-1/ctmc-strict-bd-1.xml