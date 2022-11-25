#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=MUMPS9
# Archivo de salida
#SBATCH --output=mumps_9.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=9
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 9 python3 mumps_9.py

