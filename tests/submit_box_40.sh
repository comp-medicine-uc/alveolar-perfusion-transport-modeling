#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=box40
# Archivo de salida
#SBATCH --output=box_40.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=40
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 40 python3 box_40.py
