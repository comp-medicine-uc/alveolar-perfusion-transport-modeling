#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=MUMPS3
# Archivo de salida
#SBATCH --output=mumps_3.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 3 python3 mumps_3.py

