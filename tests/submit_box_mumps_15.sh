#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=boxm15
# Archivo de salida
#SBATCH --output=box_mumps_15.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 15 python3 box_mumps_15.py
