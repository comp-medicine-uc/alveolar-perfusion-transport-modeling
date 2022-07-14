#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=nofsub_sim
# Archivo de salida
#SBATCH --output=noflatsub.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=20
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 20 python3 no_flat_sub_main.py

