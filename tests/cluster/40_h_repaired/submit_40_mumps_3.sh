#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=40mu3
# Archivo de salida
#SBATCH --output=40_mumps_3.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 3 python3 40_mumps_3.py

