#!/bin/bash

# Nombre del trabajo
#SBATCH --job-name=BoxMesh
# Archivo de salida
#SBATCH --output=box_mesh_test_1_job.txt
# Partición (Cola de y trabajo)
#SBATCH --partition=512x1024
# Solicitud de cpus
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=biherrera@uc.cl
#SBATCH --mail-type=ALL

mpirun -n 1 python3 box_mesh_test.py

