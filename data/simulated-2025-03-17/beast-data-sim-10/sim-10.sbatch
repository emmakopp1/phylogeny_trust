#!/bin/sh
# Fichier sim-10.SBATCH
#SBATCH --nodes=1
#SBATCH -c 20
#SBATCH --output=sim-10.out

~/beast/bin/beast -overwrite /home/users/kopp/work/simulated-2025-03-17/beast-data-sim-10/ctmc-strict-bd-10.xml