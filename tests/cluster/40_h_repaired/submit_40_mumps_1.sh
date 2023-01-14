#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=40mu1
# Archivo de salida
#SBATCH --output=40_mumps_1_5cpu.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 5 python3 40_mumps_1.py

