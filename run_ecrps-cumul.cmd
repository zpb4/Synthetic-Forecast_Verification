#!/bin/bash

#SBATCH -t 100:00:00
#SBATCH --job-name=ecrps-c
#SBATCH -p normal
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --nodelist=c0014
#SBATCH --output=ecrps-c.txt
#SBATCH --ntasks-per-node=1

source ~/py-env/bin/activate

python3 ./src/calc-ecrps-cumul.py