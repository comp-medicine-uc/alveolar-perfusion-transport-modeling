__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os
import time
sys.path.append(os.getcwd()[:-6])
print("Relative path: ", os.getcwd()[:-6])
import dolfin
from fenics import *
from src.model import PerfusionGasExchangeModel
from src.params import params

def run_model(name, fname, solver, preconditioner, boxmesh, max_dims, min_dims, n_jobs=1):
    print(f"Running model on {n_jobs} jobs.")
    print("Imported src files")
    print("Starting...")
    folder = fname + "/" + name
    path = os.path.join("../results-data", folder)
    model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='gmres', f_dim = 2, vecf_dim=1)
    print("Model initialised")
    
    model.import_mesh(
        os.path.join("../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
        periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
    )
    print("Mesh imported")
    
    model.mesh = boxmesh
    model.mesh = dolfin.refine(model.mesh)
    print("Mesh refined")
    
    print("Starting (P) simulation")
    model.sim_p(save=True, meshtype="tkd")
    print("(P) simulation done")
    
    print("Starting (T) simulation")
    x = model.sim_t(hb=False, save=True, solver="bicgstab")
    print("Finished (linear) guess generation")
    solution = model.sim_t(hb=True, save=True, guess=x, solver=solver, preconditioner=preconditioner)
    print("Done")


amount = 3
side_length = 5

print(f"Comienza la iteración con amount = {str(amount)} y side_length = {str(side_length)}.")

max_dims = [side_length, side_length, side_length]
min_dims = [0,0,0]

boxmesh = BoxMesh(Point(0,0,0), Point(side_length,side_length,side_length), amount, amount, amount)

name = "amount_" + str(amount)
fname = "hypre_edge_" + str(side_length)

run_model(name, fname, "bicgstab", "hypre_amg", boxmesh, max_dims, min_dims)