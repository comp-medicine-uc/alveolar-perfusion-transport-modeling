'''Main file for alveolar perfusion and gas exchange simulation in 
a flattened rat lung RVE mesh.
'''

__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os

# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-26])
print("Relative path: ", os.getcwd()[:-26])

import dolfin
from src.model import PerfusionGasExchangeModel
from src.params import params
    
print("Imported src files")
print("Starting...")
folder = "meshio_np"
path = os.path.join("../../../results-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='gmres')
print("Model initialised")
model.import_mesh(
    os.path.join("../../../raw-data/rve_40_flat/meshio", "rve_40_flat.xdmf"), type="xdmf", 
    periodic=False
)
print("Mesh imported")
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
# print("Starting (T) simulation")
# x = model.sim_t(hb=False, save=True)
# print("Finished (linear) guess generation")
# solution = model.sim_t(hb=True, save=True, guess=x)
# print("Done")