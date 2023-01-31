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


class Inlet(dolfin.SubDomain):
    def inside(self, x, on_boundary):
        return (x[0]<5-DOLFIN_EPS or x[0]>40-5+DOLFIN_EPS)

print("Imported src files")
print("Starting...")

folder = "in_prefined"
path = os.path.join("../../../results-data", folder)

model = PerfusionGasExchangeModel(folder_path=path, params=params, solver='gmres', f_dim = 2, vecf_dim=1)

max_dims = [39.894, 39.896, 39.892]
min_dims = [0.099, 0.095, 0.105]

print("Model initialised")

model.import_mesh(
    os.path.join("../../../raw-data/40_h_repaired", "40_h_repaired.xdmf"), type="xdmf", 
    periodic=False, max_dims=max_dims, min_dims=min_dims, tol=0.1
)

print("Mesh imported")

cell_markers = MeshFunction("bool", model.mesh, model.mesh.topology().dim())
cell_markers.set_all(False)

print("Defined cell markers")

inlet = Inlet()
inlet.mark(cell_markers, True)

# cell_marked = File(model.folder_path+'/bnd/cell_markers.pvd')
# cell_marked << cell_markers

new_mesh = dolfin.refine(model.mesh, cell_markers)
new_mesh_file = File(model.folder_path+'/bnd/p_refined.pvd')
new_mesh_file << new_mesh

print("Finished partial mesh refination")

model.mesh = new_mesh
print("Mesh refined")
print("Starting (P) simulation")
model.sim_p(save=True, meshtype="tkd")
print("(P) simulation done")
print("Starting (T) simulation")
x = model.sim_t(hb=False, save=True, solver="bicgstab")
print("Finished (linear) guess generation")
solution = model.sim_t(hb=True, save=True, guess=x, solver="bicgstab", preconditioner="default")
print("Done")