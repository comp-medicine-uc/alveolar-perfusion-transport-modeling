#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=fBigBoxMesh
# Archivo de salida
#SBATCH --output=fine_50_box_mesh.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=15
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 15 python3 big_box_mesh.py

