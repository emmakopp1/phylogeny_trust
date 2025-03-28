#!/bin/sh
# Fichier sim-13.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-13.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-13/ctmc-strict-bd-13.xml