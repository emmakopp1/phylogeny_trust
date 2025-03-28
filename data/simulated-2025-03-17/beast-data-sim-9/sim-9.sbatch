#!/bin/sh
# Fichier sim-9.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-9.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-9/ctmc-strict-bd-9.xml