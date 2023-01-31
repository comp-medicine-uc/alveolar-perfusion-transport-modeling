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

slowfolder = "40_h_repaired_slow"
slowpath = os.path.join("../../../results-data", slowfolder)
slowmodel = PerfusionGasExchangeModel(folder_path=slowpath, params=params)

slowmodel.params['u_in'] = 0.2
print(f"Slow model inlet velocity = {str(slowmodel.params['u_in'])}")

max_dims = [39.894, 39.896, 39.892]
min_dims = [0.099, 0.095, 0.105]

print("Slow model initialised")
slowmodel.import_mesh(
    os.path.join("../../../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)
print("Slow - mesh imported")

print("Starting slow (P) simulation")
slowmodel.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")

print("Starting slow (T) simulation")
print(f"Slow u_in value = {str(slowmodel.params['u_in'])}")
slowx = slowmodel.sim_t(hb=False, save=True)
print("Finished (linear) slow guess generation")
print("Started nonlinear slow solution")
slowsolution = slowmodel.sim_t(hb=True, save=True, guess=slowx, solver="bicgstab", preconditioner="default")
print("Finished nonlinear slow solution")

 
folder = "40_h_repaired"
path = os.path.join("../../../results-data", folder)
fastmodel = PerfusionGasExchangeModel(folder_path=path, params=params)

print("Fast model initialised")
fastmodel.import_mesh(
    os.path.join("../../../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)
      
fastmodel.params['u_in'] = 200
print("Fast - mesh imported")
print(f"Fast u_in value = {str(fastmodel.params['u_in'])}")
print("Starting fast (P) simulation")
fastmodel.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")

print("Starting fast (T) simulation without linear guess, with guess from slow solution.")
# print(f"Old u_in value = {str(model.params['u_in'])}")
# x = model.sim_t(hb=False, save=True)
# print("Finished (linear) guess generation")
solution = fastmodel.sim_t(hb=True, save=True, guess=slowsolution)
print("Done")