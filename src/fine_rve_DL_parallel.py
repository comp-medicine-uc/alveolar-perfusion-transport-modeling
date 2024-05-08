import numpy as np
import matplotlib.pyplot as plt
import pyvista
import ufl
import time
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx
from dolfinx import fem, mesh, plot, nls, log, io
from dolfinx import cpp as _cpp
import meshio
import os
import sys
import csv

sys.path.append('../')

from src.utils import *
from src.model import *
from src.diffusing_capacity_utils import *

'''
Experiments: control / emph1 / emph2 / emph3
'''

# Experiment info

control_info = [
    '../raw-data/fine_meshes/fine_rve_225_control.xdmf',
    "../results-data/fine_scale",
    "fine_rve_225_control"
]

emph_1_info = [
    '../raw-data/fine_meshes/fine_rve_225_1.xdmf',
    "../results-data/fine_scale",
    "fine_rve_225_1"
]

emph_2_info = [
    '../raw-data/fine_meshes/fine_rve_225_2.xdmf',
    "../results-data/fine_scale",
    "fine_rve_225_2"
]

emph_3_info = [
    '../raw-data/fine_meshes/fine_rve_225_3.xdmf' ,
    "../results-data/fine_scale",
    "fine_rve_225_3"
]

porosities = [0.621, 0.711, 0.794, 0.867]

all_info = [control_info, emph_1_info, emph_2_info, emph_3_info]

# Experiment number (locally, I can only run one case at a time)
i = int(sys.argv[1])

# Simulation

if MPI.COMM_WORLD.Get_rank() == 0:
    print("------------------------------------")
    print(f"Test: {all_info[i][2]}")
    print("------------------------------------")

O2_absorbed_value, CO2_released_value, DL_O2_value, DL_CO2_value, porosity = rve_test(all_info[i][0],
                                                                                        all_info[i][1],
                                                                                        all_info[i][2],
                                                                                        N=25,
                                                                                        verbose=False)

if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"################ Results for i = {i} ################")
    print(f"Diffusing capacity for oxygen DL_O2: {DL_O2_value} mL/min/mmHg")
    print(f"Diffusing capacity for carbon dioxide DL_CO2: {DL_CO2_value} mL/min/mmHg")

# Saving DLs to csv files

with open("./csv-results/"+all_info[i][2]+"_DL_O2.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow([porosity])
    writer.writerow([DL_O2_value])

with open("./csv-results/"+all_info[i][2]+"_DL_CO2.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow([porosity])
    writer.writerow([DL_CO2_value])

## Backup - mesh raw info

# mesh225_normal
# mesh_path = '../raw-data/mesh225_normal/rve_225_normal.xdmf'
# results_folder = r"../results-data/"
# exp_folder = "rve_225_normal"
# max_dims = [226.4988, 225.1387, 225.2962]
# min_dims = [-1.4082, 0.0866, 0.03899]

# mesh225_dilated
# mesh_path = '../raw-data/mesh225_dilated/rve_225_dilated.xdmf'
# results_folder = r"../results-data/"
# exp_folder = "rve_225_dilated"
# max_dims = [226.3880, 224.8622, 225.0170]
# min_dims = [-1.14634, -0.0851, 0.2545]

# mesh225_eroded
# mesh_path = '../raw-data/mesh225_eroded/rve_225_eroded.xdmf'
# results_folder = r"../results-data/"
# exp_folder = "rve_225_eroded"
# max_dims = [226.01744, 224.76256, 224.7918]
# min_dims = [-1.2616, 0.4544, 0.2872]

# mesh225_emphysema_control (porosity = 0.6211)
# mesh_path = '../raw-data/mesh225_emphysema_control/rve_225_emphysema_control.xdmf'
# results_folder = r"../results-data/"
# exp_folder = "rve_225_emphysema_control"
# max_dims = [227.47181425, 224.96266350000002, 225.09152999999998]
# min_dims = [-2.40614125, 0.0489555, -0.058911750000000006]

# mesh225_emphysema_1 (porosity = 0.7107)
# mesh_path = '../raw-data/mesh225_emphysema_1/rve_225_emphysema_1.xdmf'
# results_folder = r"../results-data/"
# exp_folder = "rve_225_emphysema_1"
# max_dims = [227.45507650000002, 225.01651275, 224.91399825]
# min_dims = [-2.5795375, 0.09627975000000001, 0.1724085]

# mesh225_emphysema_2 (porosity = 0.7942)
# mesh_path = '../raw-data/mesh225_emphysema_2/rve_225_emphysema_2.xdmf'
# results_folder = r"../results-data/"
# exp_folder = "rve_225_emphysema_2"
# max_dims = [227.54818725585938, 224.87452697753906, 224.80992126464844]
# min_dims = [-2.392364978790283, -0.0034179999493062496, 0.027458999305963516]

# mesh225_emphysema_3 (porosity = 0.8660)
# mesh_path = '../raw-data/mesh225_emphysema_3/rve_225_emphysema_3.xdmf' 
# results_folder = r"../results-data/"
# exp_folder = "rve_225_emphysema_3"
# max_dims = [227.3523712158203, 224.85842895507812, 224.78347778320312]
# min_dims = [-2.3110430240631104, 0.07416199892759323, -0.025009000673890114]