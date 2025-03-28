#!/bin/sh
# Fichier sim-15.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-15.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-15/ctmc-strict-bd-15.xml