#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=MUMPS6
# Archivo de salida
#SBATCH --output=mumps_6.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 6 python3 mumps_6.py

