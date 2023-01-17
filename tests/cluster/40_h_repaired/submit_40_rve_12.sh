#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=40rve12_1task
# Archivo de salida
#SBATCH --output=40_rve_12_1task.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 1 python3 40_rve_12.py

