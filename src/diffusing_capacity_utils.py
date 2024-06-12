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
sys.path.append('./')

from src.utils import *
from src.model import *

### Physical properties

def O2_absorbed(model, domain, ds, air_tag = 3):

    factor = 22.4 * 60 # mL/mmol * s/min
    rve_volume = 225**3 # um3
    lung_volume = 2500E12 # um3
    volume_factor = lung_volume / rve_volume
    print(f"volume factor = {volume_factor}")

    f = model.t_params['d_ba_O2'] * model.dash_params['beta_O2'] * (1/model.t_params['h_ba']) * (model.t_params['p_O2_air'] - model.p_O2)

    integral = calculate_integral_over_surface(domain, ds, [air_tag], f)
    
    return factor * volume_factor * integral

def DL_O2(model, O2_absorbed):

    return O2_absorbed / np.abs(model.t_params['p_O2_air'] - model.t_params['p_O2_in'])

def CO2_released(model, domain, ds, air_tag = 3):
    
    factor = 22.4 * 60 # mL/mmol * s/min
    rve_volume = 225**3 # um3
    lung_volume = 2500E12 # um3
    volume_factor = lung_volume / rve_volume
    print(f"volume factor = {volume_factor}")

    f = -1 * model.t_params['d_ba_CO2'] * model.dash_params['beta_CO2'] * (1/model.t_params['h_ba']) * (model.t_params['p_CO2_air'] - model.p_CO2)

    integral = calculate_integral_over_surface(domain, ds, [air_tag], f)
    
    return factor * volume_factor * integral

def DL_CO2(model, CO2_released):

    return CO2_released / np.abs(model.t_params['p_CO2_air'] - model.t_params['p_CO2_in'])

def rve_test(mesh_path, results_folder, exp_folder, max_dims=[0,0,0], min_dims=[0,0,0], uin=250, N=8, verbose=True):

    with io.XDMFFile(MPI.COMM_WORLD, mesh_path, 'r') as xdmf:
        domain = xdmf.read_mesh(name='Grid')

    model = PerfusionGasExchangeModel(mesh_path=mesh_path,
                                    results_path=results_folder,
                                    exp_path=exp_folder,
                                    params=None,
                                    ksp_type='cg', pc_type='lu', pc_factor_mat_solver_type='mumps')

    model.Setup(domain, atol = 0.1, imported=True)
    model.parameter_setup()

    model.p_params["uin"] = uin

    p, u = model.Perfusion(domain, plot=False, save=True, verbose=verbose)

    p_val = 0
    guess = model.GasExchange(domain, guess=None, save=False, 
                            plot=False, p_val = p_val, 
                            postprocess=False, plot_lines=False, verbose=verbose)
    
    for i in range(N):
        if i != N-1:
            p_val += 1/N
            guess = model.GasExchange(domain, guess=guess, save=False, 
                                    plot=False, p_val = p_val, 
                                    postprocess=False, plot_lines=False, verbose=verbose)
        else:
            p_val = 1
            guess = model.GasExchange(domain, guess=guess, save=True, 
                                    plot=False, p_val = p_val, 
                                    postprocess=True, plot_lines=False, verbose=verbose)
            
    ## Output measures

    O2_absorbed_value = O2_absorbed(model, domain, model.ds, 3)
    CO2_released_value = CO2_released(model, domain, model.ds, 3)
    DL_O2_value = DL_O2(model, O2_absorbed_value)
    DL_CO2_value = DL_CO2(model, CO2_released_value)

    return O2_absorbed_value, CO2_released_value, DL_O2_value, DL_CO2_value, model.porosity