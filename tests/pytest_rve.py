
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

# mesh_path = '../raw-data/mesh225/rve_225um.xdmf'
mesh_path = '../raw-data/mesh225_dilated/rve_225_dilated.xdmf'

with io.XDMFFile(MPI.COMM_WORLD, mesh_path, 'r') as xdmf:
    domain = xdmf.read_mesh(name='Grid')

results_folder = r"../results-data/"
exp_folder = "rve_225_dilated"

model = PerfusionGasExchangeModel(mesh_path=mesh_path,
                                  results_path=results_folder,
                                  exp_path=exp_folder,
                                  params=None)

# max_dims = [76.07255, 74.9075,  74.99651]
# min_dims = [-1.12887,  0.07919,  0.09616]

# max_dims = [226.4988, 225.138, 225.2962]
# min_dims = [-1.4086, 0.08677, 0.0389]

# mesh225_dilated
max_dims = [226.3880, 224.8622, 225.0170]
min_dims = [-1.14634, -0.0851, 0.2545]

# mesh225_normal
# max_dims = [226.4988, 225.1387, 225.2962]
# min_dims = [-1.4082, 0.0866, 0.03899]

# mesh225_eroded
# max_dims = [226.01744, 224.76256, 224.7918]
# min_dims = [-1.2616, 0.4544, 0.2872]

model.Setup(domain, atol = 1E-3, max_dims=max_dims, min_dims=min_dims, imported=True)
model.parameter_setup()

model.p_params["uin"] = 250

p, u = model.Perfusion(domain, plot=False, save=True)

p_val = 0
guess = model.GasExchange(domain, guess=None, save=False, 
                          plot=False, p_val = p_val, 
                          postprocess=False, plot_lines=False)
N = 10
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