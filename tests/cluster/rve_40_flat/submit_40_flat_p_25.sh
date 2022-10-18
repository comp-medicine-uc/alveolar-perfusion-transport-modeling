#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=flat40p
# Archivo de salida
#SBATCH --output=flat40p25.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=25
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 25 python3 40_flat_main_p_25.py

