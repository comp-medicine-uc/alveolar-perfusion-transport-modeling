#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=BigBoxMesh
# Archivo de salida
#SBATCH --output=big_box_mesh.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 4 python3 big_box_mesh.py

