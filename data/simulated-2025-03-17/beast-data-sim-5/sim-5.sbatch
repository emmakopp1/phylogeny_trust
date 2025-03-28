#!/bin/sh
# Fichier sim-5.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-5.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-5/ctmc-strict-bd-5.xml