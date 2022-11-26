#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=BiCGSt_2
# Archivo de salida
#SBATCH --output=bicgstab_2.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 2 python3 bicgstab_2.py

