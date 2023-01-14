#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=multiBox
# Archivo de salida
#SBATCH --output=multiBox.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -c 10 python3 multi.py
