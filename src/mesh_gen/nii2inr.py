# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:21:07 2022

@author: angus
"""

import os, sys
import numpy as np
import nibabel as nib

#################
### INR IMAGE ###
#################

# Writing an INR file is necessary in order to communicate the image into CGAL
# note that the header is written conventionally and that the data itself is 
# written as a bytes.

def writeINR(inrfile, field, dx, dy, dz, Type="float32", CPU="decm"):
	
	s = len(field.shape)//4  # 1 vfield 0 sfield 
	
	if s==1:
		vdim=3
	elif s==0:
		vdim=1  

	header = "#INRIMAGE-4#{\n" 
	header +="XDIM="+str(field.shape[0+s])+"\n"  # x dimension 
	header +="YDIM="+str(field.shape[1+s])+"\n"  # y dimension 
	header +="ZDIM="+str(field.shape[2+s])+"\n"  # z dimension 
	header +="VDIM="+str(vdim)+"\n"  
	header +="VX="+str(dx)+"\n"  # voxel size in x 
	header +="VY="+str(dy)+"\n"  # voxel size in y 
	header +="VZ="+str(dz)+"\n"  # voxel size in z 
	Scale = -1
	Type_np = Type.lower()

	if(Type_np in ['float32','single','float','float_vec']):
		Type = "float"
		Pixsize = "32 bits"
	elif(Type_np in ['float64','double','double_vec']):
		Type = "float"
		Pixsize = "64 bits"
	elif(Type_np in ['uint8']):
		Type = "unsigned fixed"
		Pixsize = "8 bits"
		Scale = "2**0"
	elif(Type_np in ['int16']):
		Type = "signed fixed"
		Pixsize = "16 bits"
		Scale = "2**0"
	elif(Type_np in ['uint16']):
		Type = "unsigned fixed"
		Pixsize = "16 bits"
		Scale = "2**0"
	else:
		print("Incorrect Data Type.")		

	header +="TYPE="+Type+"\n"  #float, signed fixed, or unsigned fixed 
	header +="PIXSIZE="+Pixsize+"\n"  #8, 16, 32, or 64 
	if Scale!=-1: header +="SCALE="+Scale+"\n"  #not used in my program 
	header +="CPU="+CPU+"\n"  
	# decm, alpha, pc, sun, sgi ; ittle endianness : decm, alpha,
	# pc; big endianness :sun, sgi
	for i in range(256-(len(header)+4)):
		header +="\n"
	header +="##}\n"
	
	file = open(inrfile, 'wb')
	file.write(header.encode('utf-8'))
	file.write(field.astype(Type_np).tobytes("F"))
	file.close()
	

print(" Transforming .nii.gz or .nii to an .inr file")

folder = sys.argv[1]
print("Target folder: %s"%folder)

segmentation_image = "rve.nii.gz"
inr_name = "rve.inr"

if not os.path.isfile(os.path.join(folder,segmentation_image)):
    print("** Segmentation image not found **")
    print("FATAL ERROR")
    exit()


# =============================================================================
# Reading the segmentation image
# =============================================================================
# Open the image in order to extract relevant dimension data that will also
# be useful for later computations.
print("Inspecting segmentation")
seg = nib.load(os.path.join(folder,segmentation_image))
print("   + Extracting data")
seg_data = np.array(seg.dataobj)
seg_data = seg_data[:,:,:,0,0]
print("Seg_data.shape = ", seg_data.shape)

sx, sy, sz = seg_data.shape
dx, dy, dz, _ = np.diagonal(seg.affine).__abs__()
print("     > Done")

# =============================================================================
# Writing an associated INR file
# =============================================================================
# INR is a poorly documentated filetype associated with CGAL. Documentation
# regarding it's structure is oddly found and it's mostly secondary sources
# (people commenting on threads). Currently it's working fine after being
# translated from Python 2.7.

# write INR
print("Writing INR")
writeINR(os.path.join(folder,inr_name), seg_data, dx, dy, dz)
print("     > Done")
