#!/bin/sh
# Fichier sim-4.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-4.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-4/ctmc-strict-bd-4.xml