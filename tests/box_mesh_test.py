__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os
import time

# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-6])
print("Relative path: ", os.getcwd()[:-6])
print("Dolfin version: ", dolfin.__version__)
import dolfin
from fenics import *
from src.model import PerfusionGasExchangeModel
from src.params import params

if not has_linear_algebra_backend("PETSc") and not has_linear_algebra_backend("Tpetra"):
    info("DOLFIN has not been configured with Trilinos or PETSc. Exiting.")
    exit()

if not has_krylov_solver_preconditioner("amg"):
    info("Sorry, this demo is only available when DOLFIN is compiled with AMG "
         "preconditioner, Hypre or ML.")
    exit()

if has_krylov_solver_method("minres"):
    krylov_method = "minres"
elif has_krylov_solver_method("tfqmr"):
    krylov_method = "tfqmr"
else:
    info("Default linear algebra backend was not compiled with MINRES or TFQMR "
         "Krylov subspace method. Terminating.")
    exit()
    
box_mesh = BoxMesh(Point(0,0,0),Point(5,5,5),10,10,10)

print("Imported src files")
print("Starting...")
folder = "box_mesh_on_cluster"
path = os.path.join("../results-data", folder)
print(path)
model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='mumps', f_dim = 2, vecf_dim=1)
print("Model initialised")
# max_dimsh = [12.005523, 12.005658, 12.005917]
# min_dimsh = [1.994504, 1.993421, 1.994056]
max_dimsh = [5,5,5]
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
    solver="mumps", preconditioner="default")
print(f"Time elapsed: {time.time()-t1}")
print("Done")