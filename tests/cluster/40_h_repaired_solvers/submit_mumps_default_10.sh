#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=MUMPSd10
# Archivo de salida
#SBATCH --output=mumps_default_10.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 10 python3 mumps_default_10.py

