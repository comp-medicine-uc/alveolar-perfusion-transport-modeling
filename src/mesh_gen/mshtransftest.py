import meshio
import os
# import pyvista

path = './'
src = 'out.mesh'
dest = 'out.vtk'

# print(os.listdir(path))

mesh = meshio.read(path+src)

print(type(mesh))
# mesh.cells o similar 
# extraer points y cells>tetra y revisar indexación


mesh.write(path + dest)


# path = r'//wsl.localhost/Ubuntu-18.04/home/bnherrerac/microCT-alveolar-meshgen/tests/'
# src = r'full_fine.h5'
# dest = r'full_fine.xdmf'

# # print(os.listdir(path))

# mesh = meshio.read(path+src)

# mesh.write(path + dest)


