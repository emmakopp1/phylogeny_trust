#!/bin/sh
# Fichier bantu.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=bantu-study.out

~/beast/bin/beast -overwrite ~/bantu/ctmct-strict-bd/ctmc-strict-bd.xml