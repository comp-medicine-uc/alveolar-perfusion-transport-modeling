#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=ICC
# Archivo de salida
#SBATCH --output=icc_40_h_repaired.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 20 python3 icc_40_h_repaired_main.py

