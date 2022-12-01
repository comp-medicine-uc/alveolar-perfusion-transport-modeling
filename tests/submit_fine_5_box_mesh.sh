#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=fBoxMesh
# Archivo de salida
#SBATCH --output=fine_5_box_mesh.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 1 python3 fine_5_box_mesh.py

