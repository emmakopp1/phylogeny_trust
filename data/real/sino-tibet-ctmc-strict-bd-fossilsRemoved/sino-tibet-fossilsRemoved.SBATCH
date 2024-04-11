#!/bin/sh
# Fichier fossilsRemoved.SBATCH
#SBATCH --nodes=1
#SBATCH -c 30
#SBATCH --output=ex.out

~/beast/bin/beast -overwrite /home/users/kopp/work/sino-tibetan/res/birth-death/sino-tibet-ctmc-strict-bd-fossilsRemoved/sino-tibet-ctmc-strict-bd-fossilsRemoved.xml