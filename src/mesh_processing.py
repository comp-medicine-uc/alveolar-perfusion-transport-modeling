import meshio
import os
import numpy as np
import sys
# import pyvista

path = '//wsl.localhost/Ubuntu-18.04/home/bnherrerac/alveolar-perfusion-transport-modeling/raw-and-results-data/rve_40_vol/'
src = 'rve_40_f.msh'

# For msh to xdmf conversion
# dest = src[:-4] + '.xdmf'
# mesh = meshio.read(path+src)
# mesh.write(path+dest)

# dest_unchanged = 'rve_40_vol_unchanged.vtu'
# dest_flat = 'rve_40_flat.vtu'

mesh = meshio.read(path+src)

#Get point data and modify it
points = mesh.points
# cells = mesh.cells
size = np.shape(points)
print(size)
print(points)
# # np.set_printoptions(threshold=sys.maxsize)
# print(cells)


print(size)
new_points = np.zeros(size)
print("new_points shape = ", np.shape(new_points))

# Get maxima and minima of node coordinates
max_x = np.max(points[:,0])
max_y = np.max(points[:,1])
max_z = np.max(points[:,2])
min_x = np.min(points[:,0])
min_y = np.min(points[:,1])
min_z = np.min(points[:,2])
print("max values =", [max_x, max_y, max_z])
print("min values =", [min_x, min_y, min_z])
side_length = max(max_x,max_y,max_z)-min(min_x,min_y,min_z)
print(side_length)

## Face flattening 

# Tolerance and displacement values
tol = 0.3
d = 1

for i in range(size[0]):
    current_point = points[i,:]
    if current_point[0]<min_x + tol:
        current_point[0] = min_x - d
    elif current_point[0]>max_x - tol:
        current_point[0] = max_x + d

    if current_point[1]<min_y + tol:
        current_point[1] = min_y - d
    elif current_point[1]>max_y - tol:
        current_point[1] = max_y + d

    if current_point[2]<min_z + tol:
        current_point[2] = min_z - d
    elif current_point[2]>max_z - tol:
        current_point[2] = max_z + d

    new_points[i,0] = current_point[0]
    new_points[i,1] = current_point[1]
    new_points[i,2] = current_point[2]

# Get node maxima for FEniCS mesh size configuration
nmax_x = np.max(new_points[:,0])
nmax_y = np.max(new_points[:,1])
nmax_z = np.max(new_points[:,2])
nmin_x = np.min(new_points[:,0])
nmin_y = np.min(new_points[:,1])
nmin_z = np.min(new_points[:,2])

print("nmax values =", [nmax_x, nmax_y, nmax_z])
print("nmin values =", [nmin_x, nmin_y, nmin_z])

# Updating mesh points
mesh.points = new_points

# # Exporting the mesh

vtu_name = 'rve_40_flat.xdmf'
mesh.write(path + vtu_name)

# # msh_name = 'rve_40_flat.msh'
# # mesh.write(path + msh_name)

# xdmf_name = 'rve_40_flat.xdmf'
# mesh.write(path + xdmf_name)

vtu_mesh = meshio.read(path + vtu_name)

off_name = 'rve_40_flat.off'
vtu_mesh.write(path + off_name)

ply_name = 'rve_40_flat.ply'
vtu_mesh.write(path + ply_name)

stl_name = 'rve_40_flat.stl'
vtu_mesh.write(path + stl_name)

# obj_name = 'rve_40_flat.obj'
# vtu_mesh.write(path + obj_name)

off_name = 'rve_40_flat.off'
vtu_mesh.write(path + off_name)

