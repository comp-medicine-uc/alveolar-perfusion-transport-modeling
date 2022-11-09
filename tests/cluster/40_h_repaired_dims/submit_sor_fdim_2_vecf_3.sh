#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=SOR
# Archivo de salida
#SBATCH --output=sor_40_h_repaired.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 20 python3 sor_fdim_2_vecf_3.py

