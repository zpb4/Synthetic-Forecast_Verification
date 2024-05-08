#!/bin/bash

#SBATCH -t 100:00:00
#SBATCH --job-name=extract
#SBATCH -p normal
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --nodelist=c0014
#SBATCH --output=extract.txt
#SBATCH --ntasks-per-node=1

source ~/py-env/bin/activate

python3 ./src/extract_syn-forecasts.py