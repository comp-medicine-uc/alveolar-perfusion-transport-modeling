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
    '../raw-data/mesh225_emphysema_control/rve_225_emphysema_control.xdmf',
    "../results-data/",
    "rve_225_emphysema_control",
    [227.47181425, 224.96266350000002, 225.09152999999998],
    [-2.40614125, 0.0489555, -0.058911750000000006]
]

emph_1_info = [
    '../raw-data/mesh225_emphysema_1/rve_225_emphysema_1.xdmf',
    "../results-data/",
    "rve_225_emphysema_1",
    [227.45507650000002, 225.01651275, 224.91399825],
    [-2.5795375, 0.09627975000000001, 0.1724085]
]

emph_2_info = [
    '../raw-data/mesh225_emphysema_2/rve_225_emphysema_2.xdmf',
    "../results-data/",
    "rve_225_emphysema_2",
    [227.54818725585938, 224.87452697753906, 224.80992126464844],
    [-2.392364978790283, -0.0034179999493062496, 0.027458999305963516]
]

emph_3_info = [
    '../raw-data/mesh225_emphysema_3/rve_225_emphysema_3.xdmf' ,
    "../results-data/",
    "rve_225_emphysema_3",
    [227.3523712158203, 224.85842895507812, 224.78347778320312],
    [-2.3110430240631104, 0.07416199892759323, -0.025009000673890114]
]


all_info = [control_info, emph_1_info, emph_2_info, emph_3_info]

# Experiment number
i = int(sys.argv[1])

# Simulation
if MPI.COMM_WORLD.Get_rank() == 0:
    print("------------------------------------")
    print(f"Test: {all_info[i][2]}")
    print("------------------------------------")

O2_absorbed_value, CO2_released_value, DL_O2_value, DL_CO2_value, porosity = rve_test(all_info[i][0], 
                                                                            all_info[i][1], 
                                                                            all_info[i][2], 
                                                                            all_info[i][3], 
                                                                            all_info[i][4], 
                                                                            N=12,
                                                                            verbose=True)

if MPI.COMM_WORLD.Get_rank() == 0:
    print(f"################ Results for i = {i} ################")
    print(f"Diffusing capacity for oxygen DL_O2: {DL_O2_value} mL/min/mmHg")
    print(f"Diffusing capacity for carbon dioxide DL_CO2: {DL_CO2_value} mL/min/mmHg")

# Saving DLs to csv files

with open("./csv-results/"+all_info[i][2]+"_DL_O2.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow([float(porosity)])
    writer.writerow([DL_O2_value])

with open("./csv-results/"+all_info[i][2]+"_DL_CO2.csv", "w") as file:
    writer = csv.writer(file)
    writer.writerow([porosity])
    writer.writerow([DL_CO2_value])
