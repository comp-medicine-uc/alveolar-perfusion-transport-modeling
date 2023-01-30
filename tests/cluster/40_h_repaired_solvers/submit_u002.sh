#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=u002
# Archivo de salida
#SBATCH --output=u_002wnl.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 5 python3 u_002.py

