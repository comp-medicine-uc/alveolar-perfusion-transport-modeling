# %%
import numpy as np
import matplotlib.pyplot as plt
import pyvista
import ufl
import time
from mpi4py import MPI
from petsc4py import PETSc
# import dolfinx
from dolfinx import fem, mesh, plot, nls, log, io
from dolfinx import cpp as _cpp
import meshio
import os
import seaborn as sns

from utils import *
from model import *

vertex = [120,8,8]
num_nodes = [100,5,5]
domain = mesh.create_box(MPI.COMM_WORLD, [[0.0, 0.0, 0.0], vertex], 
                         num_nodes, mesh.CellType.tetrahedron) 

results_folder = r"../../results-data/"
exp_folder = "slab_test_dim1"

model = PerfusionGasExchangeModel(mesh_path=None,
                                  results_path=results_folder,
                                  exp_path=exp_folder,
                                  params=None)

model.Setup(domain, max_dims=vertex, imported=False)
model.parameter_setup()

p, u = model.Perfusion(domain, plot=False, save=True)

p_val = 0
guess = model.GasExchange(domain, guess=None, save=False, 
                          plot=False, p_val = p_val, 
                          postprocess=False, plot_lines=False)
N = 10

t1 = time.time()
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
        
print(f"\n Problem solved using element type {model.element_type}, solver {model.ksp_type} and preconditioner {model.pc_type}")
print(f"Time elapsed: {time.time()-t1}")
print(f"Ramping amount = {N}")