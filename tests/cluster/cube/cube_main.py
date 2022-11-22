'''Main file for alveolar perfusion and gas exchange simulation in 
a flattened rat lung RVE mesh.
'''

__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os

# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-19])
print("Relative path: ", os.getcwd()[:-19])

import dolfin
from src.model import PerfusionGasExchangeModel
from src.params import params
    
print("Imported src files")
print("Starting...")
folder = "cube"
path = os.path.join("../../../results-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)

max_dims = [52.023392, 52.024105, 52.023331]
min_dims = [1.9771, 1.977327, 1.976733]

print("Model initialised")
model.import_mesh(
    os.path.join("../../../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)
print("Mesh imported")
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
print("Starting (T) simulation")
x = model.sim_t(hb=False, save=True)
print("Finished (linear) guess generation")
solution = model.sim_t(hb=True, save=True, guess=x, preconditioner='sor')
print("Done")