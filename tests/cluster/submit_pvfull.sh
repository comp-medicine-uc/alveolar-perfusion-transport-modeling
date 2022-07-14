#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=pvolflat_sim
# Archivo de salida
#SBATCH --output=pvolflat.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 20 python3 x_flattened_full_main.py

