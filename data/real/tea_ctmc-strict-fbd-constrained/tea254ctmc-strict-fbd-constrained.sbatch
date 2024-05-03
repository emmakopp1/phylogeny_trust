#!/bin/sh
# Fichier tea254ctmc-strict-fbd-constrained
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=tea254ctmc-strict-fbd-constrained.out

~/beast/bin/beast -overwrite /home/users/kopp/work/tea_language_article/tea254ctmc-strict-fbd-constrained/tea254ctmc-strict-fbd-constrained.xml