#!/bin/sh
# Fichier sim-12.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-12.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-12/ctmc-strict-bd-12.xml