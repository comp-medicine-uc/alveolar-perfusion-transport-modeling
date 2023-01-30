#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=u2
# Archivo de salida
#SBATCH --output=u_2.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 5 python3 u_2.py

