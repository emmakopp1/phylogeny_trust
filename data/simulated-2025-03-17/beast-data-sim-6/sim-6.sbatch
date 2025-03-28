#!/bin/sh
# Fichier sim-6.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-6.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-6/ctmc-strict-bd-6.xml