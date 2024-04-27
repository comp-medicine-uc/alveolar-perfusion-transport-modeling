import numpy as np
import trimesh
import tetgen
import scipy.io as sio
import pyvista as pv
import meshio
import pymeshfix as pfix

# Triangular mesh load
# Nodes and triangular elements created using iso2mesh in MATLAB, 
# exported via .mat files

nodes = sio.loadmat('../raw-data/emphysema_nodes_elems/fine_nodes_225_control.mat')
# print(nodes)
nodes = nodes['mod_nodes']
elems = sio.loadmat('../raw-data/emphysema_nodes_elems/fine_elems_225_control.mat')
elems = elems['mod_elems']
elems -= 1 # MATLAB to Python indexing conversion
print(nodes.shape)
print(elems.shape)

# Mesh cleaning (removes degeneracies, self-intersections and isolated nodes)
clean_nodes, clean_elems = pfix.clean_from_arrays(nodes, elems)

# Checks whether original mesh was already clean
print(np.array_equal(nodes, clean_nodes))
print(np.array_equal(elems, clean_elems))

# Trimesh object generation, not used so far
mesh = trimesh.Trimesh(vertices=clean_nodes,
      faces = clean_elems)

# TetGen object, for tetrahedralization
tetra_mesh = tetgen.TetGen(clean_nodes, clean_elems)

# Unused
# tetra_mesh = tetra_mesh.make_manifold()

# TetGen tetrahedralization
tnodes, telems = tetra_mesh.tetrahedralize(order=1, minratio=1.1, quality=True)

# Plots 3D mesh using PyVista
# tetra_mesh.grid.plot(show_edges=True)

# Returns max and min values for all directions, to copy them to FEniCS source code
max_x = np.max(tnodes[:,0])
max_y = np.max(tnodes[:,1])
max_z = np.max(tnodes[:,2])
min_x = np.min(tnodes[:,0])
min_y = np.min(tnodes[:,1])
min_z = np.min(tnodes[:,2])
print("max values =", [max_x, max_y, max_z])
print("min values =", [min_x, min_y, min_z])

# Writes vtu file, reads it, and generates .xdmf file using meshio

outname = 'fine_rve_225_' + 'control'

tetra_mesh.write('../raw-data/' + outname + '/' + outname + '.vtu')
# mesh = meshio.read('../raw-data/' + outname + '/' + outname + '.vtu')
# mesh.write('../raw-data/' + outname + '/' + outname + '.xdmf')

# print("tnodes=\n", tnodes)
# print("telems=\n", telems)
# print("nodes=\n", nodes)

# print("tnodes shape=", np.shape(tnodes))
# print("telems shape=", np.shape(telems))
# print("nodes shape=", np.shape(nodes))