#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=40tfq1
# Archivo de salida
#SBATCH --output=40_tfqmr_def_1.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 1 python3 40_tfqmr_def_1.py

