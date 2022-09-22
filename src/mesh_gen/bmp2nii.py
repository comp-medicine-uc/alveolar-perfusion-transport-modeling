# -*- coding: utf-8 -*-

import SimpleITK as sitk
import numpy as np
from PIL import Image


# tuple_of_arrs = ()
list_of_imgs = []
for i in range(0,113):
    current_name = "../../raw-and-results-data/imgs/2018-10-12 pulmon _4.5um_rec0000_voi_0" + str(645+i) + ".bmp"
    current_img = sitk.ReadImage(current_name, imageIO="BMPImageIO")
    # current_array = (np.array(current_img)>20)*1
    # tuple_of_arrs = tuple_of_arrs + (current_array,)
    list_of_imgs.append(current_img)
    out_name = "img"+str(i)+".nii.gz"

series = sitk.JoinSeries(list_of_imgs)
sitk.WriteImage(series, "all_img.nii.gz")

# all_seg = np.dstack(tuple_of_arrs)
# print(all_seg[:10,:10,0])

# imagen = sitk.GetImageFromArray(all_seg)

# for img in tuple_of_imgs:
#     sitk.WriteImage(img, 'intento.nii', True)





