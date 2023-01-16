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
folder = "rve_40_12"
path = os.path.join("../../../results-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='gmres', f_dim = 2, vecf_dim=1)

max_dims = [39.9592, 40.0040, 39.9904]
min_dims = [0.02681, -0.01336, 0.03181]

print("Model initialised")
model.import_mesh(
    os.path.join("../../../raw-data/40_rve_12", "rve_40_12.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)
print("Mesh imported")
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
print("Starting (T) simulation")
x = model.sim_t(hb=False, save=True, solver="bicgstab")
print("Finished (linear) guess generation")
solution = model.sim_t(hb=True, save=True, guess=x, solver="mumps", preconditioner=None)
print("Done")