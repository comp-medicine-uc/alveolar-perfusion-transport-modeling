import numpy as np
import matplotlib.pyplot as plt
import trimesh
import tetgen

original_shape = (50,50)
pad = 5

square = np.ones(original_shape)
padded_square = np.pad(square, pad, mode='constant', constant_values=0)

cube = np.zeros((original_shape[0]+2*pad, original_shape[1]+2*pad, original_shape[1]+2*pad))
for i in range(0,4):
    cube[:,:,i] = np.zeros(np.shape(cube[:,:,i]))

for i in range(original_shape[1]+2*pad-pad, original_shape[1]+2*pad):
    cube[:,:,i] = np.zeros(np.shape(cube[:,:,i]))

for i in range(5, original_shape[1]+2*pad-pad-1):
    cube[:,:,i] = padded_square


# fig, ax = plt.subplots()
# ax.imshow(padded_square)
# plt.show()