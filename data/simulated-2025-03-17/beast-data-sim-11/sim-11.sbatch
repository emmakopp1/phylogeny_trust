#!/bin/sh
# Fichier sim-11.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-11.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-11/ctmc-strict-bd-11.xml