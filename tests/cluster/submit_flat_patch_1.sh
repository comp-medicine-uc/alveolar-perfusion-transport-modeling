#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=flatp1
# Archivo de salida
#SBATCH --output=flatp1.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 20 python3 flat_patch_main_1.py

