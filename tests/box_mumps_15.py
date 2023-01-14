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
print("Imported src files")

def run_model(name, fname, solver, preconditioner, max_dims, min_dims, side_length, amount):
    print("Starting...")
    folder = fname + "/" + name
    path = os.path.join("../results-data", folder)
    model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='gmres', f_dim = 2, vecf_dim=1)
    print("Model initialised")
    
    model.import_mesh(max_dims=max_dims, min_dims=min_dims, tol=0.1, box_side_length = side_length, box_nodes = amount)
    print("Mesh initialized")
    
    model.mesh = dolfin.refine(model.mesh) # Esto duplica la cantidad de nodos
    print("Mesh refined")
    
    print("Starting (P) simulation")
    model.sim_p(save=True, meshtype="tkd")
    print("(P) simulation done")
    
    print("Starting linear (T) simulation")
    x = model.sim_t(hb=False, save=True, solver="bicgstab")
    print("Finished (linear) guess generation")
    
    print("Starting nonlinear (T) simulation")
    solution = model.sim_t(hb=True, save=True, guess=x, solver=solver, preconditioner=preconditioner)
    print("Generated non-linear solution")
    print("Done")


####################### VARIAR ESTOS DOS PARÁMETROS

amount = 50
side_length = 25

####################### Ejecución

print(f"Comienza la iteración con amount = {str(amount)} y side_length = {str(side_length)}.")

max_dims = [side_length, side_length, side_length]
min_dims = [0,0,0]

name = "amount_" + str(amount)
fname = "box_gmres_15_edge_" + str(side_length)

run_model(name, fname, "gmres", "default", max_dims, min_dims, side_length, amount)