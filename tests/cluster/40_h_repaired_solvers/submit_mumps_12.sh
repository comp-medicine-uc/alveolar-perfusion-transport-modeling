#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=MUMPS12
# Archivo de salida
#SBATCH --output=mumps_12.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 12 python3 mumps_12.py

