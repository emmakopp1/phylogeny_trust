#!/bin/sh
# Fichier sim-2.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-2.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-2/ctmc-strict-bd-2.xml