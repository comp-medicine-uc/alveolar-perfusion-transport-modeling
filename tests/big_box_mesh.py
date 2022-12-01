__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os
import time
# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-6])
print("Relative path: ", os.getcwd()[:-6])

import dolfin
# print("Dolfin version: ", dolfin.__version__)
from fenics import *
from src.model import PerfusionGasExchangeModel
from src.params import params
    
box_mesh = BoxMesh(Point(0,0,0),Point(50,50,50),100,100,100)

# 33^3 ~ 36000 cells, which is approximately the same as 40_h_repaired data

print("Imported src files")
print("Starting...")
folder = "big_box_mesh"
path = os.path.join("../results-data", folder)
print(path)
model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='mumps', f_dim = 2, vecf_dim=1)
print("Model initialised")
max_dimsh = [50,50,50]
min_dimsh = [0,0,0]

model.import_mesh(
    os.path.join("../raw-data/cube", "cube.xdmf"), type="xdmf", 
    periodic=False, max_dims = max_dimsh, min_dims = min_dimsh, tol=0.1
)
print("Mesh imported")
# model.mesh = dolfin.refine(model.mesh)
model.mesh = box_mesh
# print("Mesh refined")
print("Starting (P) simulation")
t_1 = time.time()
model.sim_p(save=True, meshtype="tkd")
print("Time elapsed: ", time.time()-t_1)
print("(P) simulation done")

print("Starting (T) simulation")
t1 = time.time()
x = model.sim_t(hb=False, save=True, solver="bicgstab")
print(f"Time elapsed: {time.time()-t1}")
print("Finished (linear) guess generation")

t1 = time.time()
solution = model.sim_t(
    hb=True, save=True, guess=x, 
    solver="bicgstab", preconditioner="default")
print(f"Time elapsed: {time.time()-t1}")
print("Done")
