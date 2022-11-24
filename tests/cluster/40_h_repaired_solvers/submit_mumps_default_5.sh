#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=MUMPSdef
# Archivo de salida
#SBATCH --output=mumps_default_1.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 5 python3 mumps_default.py

