#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=opt_3
# Archivo de salida
#SBATCH --output=opt_3.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 10 python3 opt_3.py

