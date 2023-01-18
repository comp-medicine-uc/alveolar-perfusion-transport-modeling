#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=40tfqa10
# Archivo de salida
#SBATCH --output=40_tfqmr_amg_10.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 10 python3 40_tfqmr_amg_10.py
