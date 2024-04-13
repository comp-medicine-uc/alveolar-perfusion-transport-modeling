
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

sys.path.append('../')

from src.utils import *
from src.model import *

def rve_test(mesh_path, results_folder, exp_folder, max_dims, min_dims, uin=250, N=8):

    with io.XDMFFile(MPI.COMM_WORLD, mesh_path, 'r') as xdmf:
        domain = xdmf.read_mesh(name='Grid')

    model = PerfusionGasExchangeModel(mesh_path=mesh_path,
                                    results_path=results_folder,
                                    exp_path=exp_folder,
                                    params=None)

    model.Setup(domain, atol = 1E-3, max_dims=max_dims, min_dims=min_dims, imported=True)
    model.parameter_setup()

    model.p_params["uin"] = uin

    p, u = model.Perfusion(domain, plot=False, save=True)

    p_val = 0
    guess = model.GasExchange(domain, guess=None, save=False, 
                            plot=False, p_val = p_val, 
                            postprocess=False, plot_lines=False)
    
    for i in range(N):
        p_val += 1/N
        if i != N-1:
            guess = model.GasExchange(domain, guess=guess, save=False, 
                                    plot=False, p_val = p_val, 
                                    postprocess=False, plot_lines=False)
        else:
            guess = model.GasExchange(domain, guess=guess, save=True, 
                                    plot=False, p_val = p_val, 
                                    postprocess=True, plot_lines=False)
            
    # ## Output measures

    O2_absorbed_value = O2_absorbed(model, domain, model.ds, 3)
    CO2_released_value = CO2_released(model, domain, model.ds, 3)
    DL_O2_value = DL_O2(model, O2_absorbed_value)
    DL_CO2_value = DL_CO2(model, CO2_released_value)

    return O2_absorbed_value, CO2_released_value, DL_O2_value, DL_CO2_value

def multi_rve_test(mesh_path_list, results_folder, exp_folder_list, max_dims_list, min_dims_list, uin=250, N=8):

    outs = []

    for mesh_path, exp_folder, i in zip(mesh_path_list, exp_folder_list, range(len(exp_folder_list))):

        O2_absorbed_value, CO2_released_value, DL_O2_value, DL_CO2_value = rve_test(mesh_path, results_folder, exp_folder, max_dims_list[i], min_dims_list[i], uin=uin, N=N)

        outs[i] = {'O2_absorbed': O2_absorbed_value,
                   'CO2_released': CO2_released_value,
                   'DL_O2': DL_O2_value,
                   'DL_CO2': DL_CO2_value}

# mesh_path = '../raw-data/mesh225/rve_225um.xdmf'
        
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
mesh_path = '../raw-data/mesh225_emphysema_3/rve_225_emphysema_3.xdmf' 
results_folder = r"../results-data/"
exp_folder = "rve_225_emphysema_3"
max_dims = [227.3523712158203, 224.85842895507812, 224.78347778320312]
min_dims = [-2.3110430240631104, 0.07416199892759323, -0.025009000673890114]

print("------------------------------------")
print(f"Test: {exp_folder}")
print("------------------------------------")

O2_absorbed_value, CO2_released_value, DL_O2_value, DL_CO2_value = rve_test(mesh_path, results_folder, exp_folder, max_dims, min_dims, N=12)

print(f"O2 absorbed: {O2_absorbed_value} mL/min")
print(f"CO2 released: {CO2_released_value} mL/min")
print(f"Diffusing capacity for oxygen DL_O2: {DL_O2_value} mL/min/mmHg")
print(f"Diffusing capacity for carbon dioxide DL_CO2: {DL_CO2_value} mL/min/mmHg")

