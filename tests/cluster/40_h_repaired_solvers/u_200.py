'''Main file for alveolar perfusion and gas exchange simulation in 
a flattened rat lung RVE mesh.
'''

__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os

# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-36])
print("Relative path: ", os.getcwd()[:-36])

import dolfin
from src.model import PerfusionGasExchangeModel
from src.params import params

print("Imported src files")
print("Starting...")
folder = "u_200"
path = os.path.join("../../../results-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='gmres', f_dim = 2, vecf_dim=1)

max_dims = [39.894, 39.895, 39.892]
min_dims = [0.099, 0.095, 0.104]

print("Model initialised")
model.import_mesh(
    os.path.join("../../../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)
print("Mesh imported")
# model.mesh = dolfin.refine(model.mesh)
# print("Mesh refined")
model.params["u_in"] = 200
print(f"Changed inlet speed u_in to {str(model.params['u_in'])}")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
print("Starting (T) simulation")
x = model.sim_t(hb=False, save=True, solver="bicgstab")
print("Finished (linear) guess generation")
solution = model.sim_t(hb=True, save=True, guess=x, solver="bicgstab", preconditioner="default")
print("Done")