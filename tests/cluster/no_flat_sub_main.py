'''Main file for alveolar perfusion and gas exchange simulations in TKD
mesh.
'''

__author__ = 'biherrera'
__email__ = 'biherrera@uc.cl'

import sys
import os
# The following line adds the directory to the path in order to cross-reference
# files in the repo
sys.path.append(os.getcwd()[:-14])
print("Relative path: ", os.getcwd()[:-14])
import dolfin
from src.model import PerfusionGasExchangeModel
from src.params import params

print("Starting...")
folder = "non-flattened-sub-RVE"
path = os.path.join("../../raw-and-results-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)
print("Model initialised")
model.import_mesh(
    os.path.join("../../raw-and-results-data/non-flattened-sub-RVE", "non_flattened_sub_RVE_fixed.xdmf"), type="xdmf", 
    periodic=False
)
print("Mesh imported")
model.mesh = dolfin.refine(model.mesh)
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
# print("Starting (T) simulation")
# x = model.sim_t(hb=False, save=False)
# solution = model.sim_t(hb=True, save=True, guess=x)
# print("Done")