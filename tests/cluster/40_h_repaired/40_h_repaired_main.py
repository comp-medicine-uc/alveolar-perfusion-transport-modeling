'''Main file for alveolar perfusion and gas exchange simulation in 
a flattened rat lung RVE mesh.
'''

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

folder = "40_h_repaired_sf"
path = os.path.join("../../../results-data", folder)
model = PerfusionGasExchangeModel(folder_path=path, params=params)

model.params['u_in'] = 0.2
print(f"Current (slow) model inlet velocity u_in = {str(model.params['u_in'])}")

max_dims = [39.894, 39.896, 39.892]
min_dims = [0.099, 0.095, 0.105]

print("Model initialised")
model.import_mesh(
    os.path.join("../../../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)
print("Mesh imported")
model.folder_path = os.path.join(path, "/slow")
      
print(f"Starting slow (P) simulation with inlet velocity u_in {str(model.params['u_in'])}")
model.sim_p(save=True, meshtype="tkd")
print("Slow (P) simulation done")
print("Starting slow (T) simulation")
print(f"Slow u_in value = {str(model.params['u_in'])}")
slowx = model.sim_t(hb=False, save=True)
print("Finished (linear) slow guess generation")
print("Started nonlinear slow solution")
slowsolution = model.sim_t(hb=True, save=True, guess=slowx, solver="bicgstab", preconditioner="default")
print("Finished nonlinear slow solution")

print("%%%%%%%%%%%%%%%%%%%%%%%%")
      
print("Starting fast simulation")
model.params['u_in'] = 200
print(f"Current (fast) model inlet velocity u_in = {str(model.params['u_in'])}")
model.folder_path = os.path.join(path, "/fast") 
      
print("Starting fast (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("Fast (P) simulation done with inlet velocity u_in = {str(model.params['u_in'])}")

print("Starting fast (T) simulation without linear guess, with guess from slow solution.")
print(f"Current (fast) model inlet velocity u_in = {str(model.params['u_in'])}")
# print(f"Old u_in value = {str(model.params['u_in'])}")
# x = model.sim_t(hb=False, save=True)
# print("Finished (linear) guess generation")
solution = model.sim_t(hb=True, save=True, guess=slowsolution)
print("Done")