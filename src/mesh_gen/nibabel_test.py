import nibabel as nb

image = nb.load("./all_img.nii.gz")
print("image shape =", image.shape)
matrix = image.affine
print("matrix shape =", matrix.shape)
for i in range(0,3):
    print([matrix[i,i]*image.shape[i]])
print(matrix)

print("loaded")