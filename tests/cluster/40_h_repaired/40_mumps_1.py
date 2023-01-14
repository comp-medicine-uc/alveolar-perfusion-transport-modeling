__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os

# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-28])
print("Relative path: ", os.getcwd()[:-28])

import dolfin
from src.model import PerfusionGasExchangeModel
from src.params import params
    
print("Imported src files")
print("Starting...")
folder = "40_gmres_1"
path = os.path.join("../../../results-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='gmres', f_dim = 2, vecf_dim=1)

max_dims = [39.894161224365234, 39.895729064941406, 39.89208984375]
min_dims = [0.09939099848270416, 0.09558500349521637, 0.1048400029540062]

print("Model initialised")
model.import_mesh(
    os.path.join("../../../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)
print("Mesh imported")
# model.mesh = dolfin.refine(model.mesh)
print("Mesh HAS NOT BEEN refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
print("Starting (T) simulation")
x = model.sim_t(hb=False, save=True, solver="bicgstab")
print("Finished (linear) guess generation")
solution = model.sim_t(hb=True, save=True, guess=x, solver="gmres", preconditioner="hypre_amg")
print("Done")