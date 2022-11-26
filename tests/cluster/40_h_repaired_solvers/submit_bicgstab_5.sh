#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=BiCGSt_5
# Archivo de salida
#SBATCH --output=bicgstab_5.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=5
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 5 python3 bicgstab_5.py

