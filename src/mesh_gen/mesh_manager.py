# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 20:51:43 2022

@author: angus
"""

import os, sys
import meshio as io
import trimesh,tetgen
import numpy as np
import nibabel as nib
from numpy.core.umath_tests import matrix_multiply 
from scipy.sparse import coo_matrix


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


def Evaluate_FlyingNodes(xyz,IEN):
    '''
    Rutina que elimina los nodos "voladores" (con masa cero). Nodos que se pueden haber generado pero que no estan conectados.
    
    Input:
        filename_npz: archivo .npz que sera evaluado.
    Output:
        filename_Output_npz: archivo .npz con nodos voladores eliminados.
    '''
    
    xyz=np.ndarray.tolist(xyz)
    IEN=np.ndarray.tolist(IEN)
    
    aux=np.zeros((len(xyz)),dtype=bool)
    aux=np.ndarray.tolist(aux)
    
    ele=len(IEN)
    
    for i in range(ele):
        for j in range(4):
            aux[IEN[i][j]]=True
    
    nodosVoladores=[]
    for i in range(len(aux)):
        if aux[i]==False:
            nodosVoladores.append(i)
    
    xyz_NEW = np.copy(xyz)
    xyz_NEW = np.ndarray.tolist(xyz_NEW)
    
    cont=0
    for j in range(len(nodosVoladores)):
        xyz_NEW=np.delete(xyz_NEW,nodosVoladores[j]-cont,0)    
        cont=cont+1
                
    IEN_NEW = np.copy(IEN)
    IEN_NEW = np.ndarray.tolist(IEN_NEW )
    
    cont1=0
    for n in range(len(nodosVoladores)):
        for i in range(ele):
            for j in range(4):
                if IEN_NEW[i][j] > nodosVoladores[n]-cont1:
                    IEN_NEW[i][j] = IEN_NEW[i][j]-1
        cont1=cont1+1
    
    IEN_NEW = np.array(IEN_NEW)
        
    return xyz_NEW, IEN_NEW

def InterpolateBSplines(filenameResults, filenameCPP, filenameRef):
    '''
    Interpolate data using cubic B-splines.
    filenameResults : .npz file containing nodal information
    filenameCPP     : file containing control point position from registration process
    filenameRef     : Image file
    
    Output: Phi (nodal deformation mapping)
    '''    
    
    #Filename Results
    npzfile = np.load(filenameResults)
    X = npzfile['xyz']
    elem = npzfile['elem']        
    
    #Reference Affine
    img = nib.load(filenameRef)
    affine=img.affine
    #Affine 0 (origin in (0,0,0) and positive spacing) 
    affine0=np.zeros((4,4))
    affine0[0:3,0:3]=np.abs(affine[0:3,0:3])
    affine0[3,3]=1.0
    affine0_inv=np.array(np.matrix(affine0)**-1)

    #Real World Coordinates (Affine0) 2 Voxel Coordinates  
    X_VoxCoord_Affine0=np.array(np.matrix(affine0_inv[0:3,0:3])*np.matrix(X).T)
    ones=np.ones((1,X_VoxCoord_Affine0.shape[1]))
    X_VoxCoord_Affine0=np.vstack((X_VoxCoord_Affine0,ones))
    
    # Voxel Coordinates 2 Real World Coordinates (Reference Image)
    X_RWC=np.matrix(affine)*np.matrix(X_VoxCoord_Affine0)
    X_RWC=np.array(X_RWC[0:3,:].T)

    #CCP and Affine of the CPP Image
    imgCPP = nib.load(filenameCPP)
    affineCPP=imgCPP.affine
    
    dataCPP=np.squeeze(np.array(imgCPP.dataobj))
    
    #CPP Affine 
    A=np.matrix(affineCPP[0:3,0:3])**-1
    b=affineCPP[0:3,3]
    
    #COORDENADAS INDICIALES CPP DE X, FLOOR Y LOCAL
    X_VoxelCoord_CPP=np.array((np.matrix(A)*np.matrix(X_RWC-b).T).T)
    X_VoxelCoord_CPP_Floor=np.floor(X_VoxelCoord_CPP)
    X_VoxelCoord_CPP_Local=X_VoxelCoord_CPP-X_VoxelCoord_CPP_Floor

    
    #FUNCIONES BSPLINE
    def B0(u):
        return (1.-u)**3./6.
    def B1(u):
        return (3.*u**3.-6.*u**2.+4.)/6.
    def B2(u):
        return (-3.*u**3.+3.*u**2.+3.*u+1.)/6.
    def B3(u):
        return u**3/6
    
    defmap=[]
    
    #GRID A CONSIDERAR EN LA INTERPOLACION EN CADA REGION (64 PUNTOS)
    i0=np.arange(0,4)   
    j0=np.arange(0,4)
    k0=np.arange(0,4)
    i0,j0,k0=np.meshgrid(i0,j0,k0)
    i0=i0.reshape((i0.shape[0]*i0.shape[1]*i0.shape[2],))
    j0=j0.reshape((j0.shape[0]*j0.shape[1]*j0.shape[2],))
    k0=k0.reshape((k0.shape[0]*k0.shape[1]*k0.shape[2],))   
        
    for p in np.arange(X_VoxelCoord_CPP.shape[0]):
       iref=X_VoxelCoord_CPP_Floor[p,0].astype(int)-1
       jref=X_VoxelCoord_CPP_Floor[p,1].astype(int)-1
       kref=X_VoxelCoord_CPP_Floor[p,2].astype(int)-1
       i=i0+iref
       j=j0+jref
       k=k0+kref   
       u=X_VoxelCoord_CPP_Local[p,0]
       v=X_VoxelCoord_CPP_Local[p,1]
       w=X_VoxelCoord_CPP_Local[p,2]
       bu=np.array([B0(u),B1(u),B2(u),B3(u)])
       bv=np.array([B0(v),B1(v),B2(v),B3(v)])
       bw=np.array([B0(w),B1(w),B2(w),B3(w)])
    
       defmap.append([sum(bu[i0]*bv[j0]*bw[k0]*dataCPP[i,j,k,0]),
                   sum(bu[i0]*bv[j0]*bw[k0]*dataCPP[i,j,k,1]),
                   sum(bu[i0]*bv[j0]*bw[k0]*dataCPP[i,j,k,2])])
       
    PHI_RWC=np.array(defmap)

    #Affine 0 (origin in (0,0,0) and positive spacing) 
    A=np.matrix(affine[0:3,0:3])**-1
    b=affine[0:3,3]
    PHI_VoxelCoord_RefImage=np.array((A*np.matrix(PHI_RWC-b).T))
    PHI=np.matrix(affine0[0:3,0:3])*np.matrix(PHI_VoxelCoord_RefImage)
    PHI=np.array(PHI.T)
    np.savez(filenameResults,xyz=X,elem=elem,phi=PHI)
    return PHI



def interpolateIntensities(npz, img_ref, img_flo):

    for path in [npz, img_ref, img_flo]:
        if not os.path.isfile(path):
            print("%s not found"%path.split("/")[-1])
            print("SHUTTING DOWN THE PROGRAM")
            exit()
    
    warpedIntensities_ee = []
    warpedIntensities_ei = []
    
    for rqst, image in [("xyz",img_ref),("phi",img_flo)]:
        
        # load image
        im = nib.load(image)
        
        # load ref image
        im_ref = nib.load(img_ref)

        # retrieve affine and data
        affine = im_ref.affine
        data = np.array(im.dataobj)
        
        # load mesh and retrieve the corresponding points
        mesh = np.load(npz)
        X = mesh[rqst]
        
        if rqst == "phi":
            # "Z-correction" TESTING 
            z_f = im.affine[2,3]
            z_r = im_ref.affine[2,3]
            dz = z_r - z_f
            X[:,2] -= dz
        
        # interpolation routine
        affine0 = np.zeros((4,4))
        affine0[0:3,0:3] = np.abs(affine[0:3,0:3]) 
        affine0[3,3] = 1.0
        affine0_inv=np.array(np.matrix(affine0)**-1)
        
        X_VoxCoord = np.array(np.matrix(affine0_inv[0:3,0:3])*np.matrix(X).T).T
        X_VoxCoord_Floor = np.floor(X_VoxCoord)
        X_VoxCoord_Local = X_VoxCoord - X_VoxCoord_Floor
        
        i0 = np.arange(0,4)
        j0 = np.arange(0,4)
        k0 = np.arange(0,4)
        i0,j0,k0 = np.meshgrid(i0,j0,k0)
        i0 = i0.reshape((i0.shape[0]*i0.shape[1]*i0.shape[2],))
        j0 = j0.reshape((j0.shape[0]*j0.shape[1]*j0.shape[2],))
        k0 = k0.reshape((k0.shape[0]*k0.shape[1]*k0.shape[2],))
        
        def B0(u):
            return (1.-u)**3./6.
        def B1(u):
            return (3.*u**3.-6.*u**2.+4.)/6.
        def B2(u):
            return (-3.*u**3.+3.*u**2.+3.*u+1.)/6.
        def B3(u):
            return u**3/6
                
        for p in np.arange(X_VoxCoord.shape[0]):
        
            iref=X_VoxCoord_Floor[p,0].astype(int)-1
            jref=X_VoxCoord_Floor[p,1].astype(int)-1
            kref=X_VoxCoord_Floor[p,2].astype(int)-1
        
            i=i0+iref
            j=j0+jref
            k=k0+kref
            
            u=X_VoxCoord_Local[p,0]
            v=X_VoxCoord_Local[p,1]
            w=X_VoxCoord_Local[p,2]
        
            bu=np.array([B0(u),B1(u),B2(u),B3(u)])
            bv=np.array([B0(v),B1(v),B2(v),B3(v)])
            bw=np.array([B0(w),B1(w),B2(w),B3(w)])
            
            if rqst == "xyz":
                warpedIntensities_ee.append(sum(bu[i0]*bv[j0]*bw[k0]*data[i,j,k]))
            else:
                warpedIntensities_ei.append(sum(bu[i0]*bv[j0]*bw[k0]*data[i,j,k]))

    return np.array(warpedIntensities_ee), np.array(warpedIntensities_ei)



def Tet4DefGradSmoothingOptiSPARSE(xyz, defmap, LM):
    """
    This subroutine computes the assemble of mass vector and force vector (matrix)
    Input:
      xyz        : array with node coordinates in reference configuration
      defmap     : array with nodal deformation mapping vector 
      LM         : connectivity matrix

    Output:
      Mlumped   : Lumped mass matrix (in vector form)
      P         : Weighted deformation-gradient vector
    """    
    me  = LM
    q   = xyz
    phi = defmap
    
    V = np.zeros((me.shape[0],4,4))
    
    qme   = q[me[:]]
    phime = phi[me[:]]
    
    ones = np.ones((me.shape[0],4))
    x = qme[:,:,0]
    y = qme[:,:,1]
    z = qme[:,:,2]
    
    T = np.column_stack((ones[:],x[:],y[:],z[:]))
    V[:,0,0]=T[:,0];V[:,0,1]=T[:,1];V[:,0,2]=T[:,2];V[:,0,3]=T[:,3];
    V[:,1,0]=T[:,4];V[:,1,1]=T[:,5];V[:,1,2]=T[:,6];V[:,1,3]=T[:,7];
    V[:,2,0]=T[:,8];V[:,2,1]=T[:,9];V[:,2,2]=T[:,10];V[:,2,3]=T[:,11]
    V[:,3,0]=T[:,12];V[:,3,1]=T[:,13];V[:,3,2]=T[:,14];V[:,3,3]=T[:,15]
    
    ve0=np.zeros(me.shape[0])
    ve0[:]=np.linalg.det(V[:])
    ve0=ve0/6.
    
    xdef=phime[:,:,0]
    ydef=phime[:,:,1]
    zdef=phime[:,:,2]
    
    T=np.column_stack((ones[:],xdef[:],ydef[:],zdef[:]))
    
    V[:,0,0]=T[:,0];V[:,0,1]=T[:,1];V[:,0,2]=T[:,2];V[:,0,3]=T[:,3];
    V[:,1,0]=T[:,4];V[:,1,1]=T[:,5];V[:,1,2]=T[:,6];V[:,1,3]=T[:,7];
    V[:,2,0]=T[:,8];V[:,2,1]=T[:,9];V[:,2,2]=T[:,10];V[:,2,3]=T[:,11]
    V[:,3,0]=T[:,12];V[:,3,1]=T[:,13];V[:,3,2]=T[:,14];V[:,3,3]=T[:,15]
    
    ve=np.zeros(me.shape[0])
    ve[:]=np.linalg.det(V[:])
    ve=ve/6.
    
    x=qme[:,:,0]
    y=qme[:,:,1]
    z=qme[:,:,2]
    
    a1=(y[:,3]-y[:,1])*(z[:,2]-z[:,1])-(y[:,2]-y[:,1])*(z[:,3]-z[:,1])      #yd[3,1]*zd[2,1]-yd[2,1]*zd[3,1]
    a2=(y[:,2]-y[:,0])*(z[:,3]-z[:,2])-(y[:,2]-y[:,3])*(z[:,0]-z[:,2])      #yd[2,0]*zd[3,2]-yd[2,3]*zd[0,2]
    a3=(y[:,1]-y[:,3])*(z[:,0]-z[:,3])-(y[:,0]-y[:,3])*(z[:,1]-z[:,3])      #yd[1,3]*zd[0,3]-yd[0,3]*zd[1,3]
    a4=(y[:,0]-y[:,2])*(z[:,1]-z[:,0])-(y[:,0]-y[:,1])*(z[:,2]-z[:,0])      #yd[0,2]*zd[1,0]-yd[0,1]*zd[2,0]
    
    b1=(x[:,2]-x[:,1])*(z[:,3]-z[:,1])-(x[:,3]-x[:,1])*(z[:,2]-z[:,1])      #xd[2,1]*zd[3,1]-xd[3,1]*zd[2,1]
    b2=(x[:,3]-x[:,2])*(z[:,2]-z[:,0])-(x[:,0]-x[:,2])*(z[:,2]-z[:,3])      #xd[3,2]*zd[2,0]-xd[0,2]*zd[2,3]
    b3=(x[:,0]-x[:,3])*(z[:,1]-z[:,3])-(x[:,1]-x[:,3])*(z[:,0]-z[:,3])      #xd[0,3]*zd[1,3]-xd[1,3]*zd[0,3]
    b4=(x[:,1]-x[:,0])*(z[:,0]-z[:,2])-(x[:,2]-x[:,0])*(z[:,0]-z[:,1])      #xd[1,0]*zd[0,2]-xd[2,0]*zd[0,1]
    
    c1=(x[:,3]-x[:,1])*(y[:,2]-y[:,1])-(x[:,2]-x[:,1])*(y[:,3]-y[:,1])      #xd[3,1]*yd[2,1]-xd[2,1]*yd[3,1]
    c2=(x[:,2]-x[:,0])*(y[:,3]-y[:,2])-(x[:,2]-x[:,3])*(y[:,0]-y[:,2])      #xd[2,0]*yd[3,2]-xd[2,3]*yd[0,2]
    c3=(x[:,1]-x[:,3])*(y[:,0]-y[:,3])-(x[:,0]-x[:,3])*(y[:,1]-y[:,3])      #xd[1,3]*yd[0,3]-xd[0,3]*yd[1,3]
    c4=(x[:,0]-x[:,2])*(y[:,1]-y[:,0])-(x[:,0]-x[:,1])*(y[:,2]-y[:,0])      #xd[0,2]*yd[1,0]-xd[0,1]*yd[2,0]
    
    
    DN=[[a1[:]/(6*ve0[:]),a2[:]/(6*ve0[:]),a3[:]/(6*ve0[:]),a4[:]/(6*ve0[:])],[b1[:]/(6*ve0[:]),b2[:]/(6*ve0[:]),b3[:]/(6*ve0[:]),b4[:]/(6*ve0[:])],[c1[:]/(6*ve0[:]),c2[:]/(6*ve0[:]),c3[:]/(6*ve0[:]),c4[:]/(6*ve0[:])]]        
    
    DN=np.array(DN).T
    
    F11=DN[:,0,0]*phime[:,0,0]+DN[:,1,0]*phime[:,1,0]+DN[:,2,0]*phime[:,2,0]+DN[:,3,0]*phime[:,3,0]
    F12=DN[:,0,1]*phime[:,0,0]+DN[:,1,1]*phime[:,1,0]+DN[:,2,1]*phime[:,2,0]+DN[:,3,1]*phime[:,3,0]
    F13=DN[:,0,2]*phime[:,0,0]+DN[:,1,2]*phime[:,1,0]+DN[:,2,2]*phime[:,2,0]+DN[:,3,2]*phime[:,3,0]
    
    F21=DN[:,0,0]*phime[:,0,1]+DN[:,1,0]*phime[:,1,1]+DN[:,2,0]*phime[:,2,1]+DN[:,3,0]*phime[:,3,1]
    F22=DN[:,0,1]*phime[:,0,1]+DN[:,1,1]*phime[:,1,1]+DN[:,2,1]*phime[:,2,1]+DN[:,3,1]*phime[:,3,1]
    F23=DN[:,0,2]*phime[:,0,1]+DN[:,1,2]*phime[:,1,1]+DN[:,2,2]*phime[:,2,1]+DN[:,3,2]*phime[:,3,1]
    
    F31=DN[:,0,0]*phime[:,0,2]+DN[:,1,0]*phime[:,1,2]+DN[:,2,0]*phime[:,2,2]+DN[:,3,0]*phime[:,3,2]
    F32=DN[:,0,1]*phime[:,0,2]+DN[:,1,1]*phime[:,1,2]+DN[:,2,1]*phime[:,2,2]+DN[:,3,1]*phime[:,3,2]
    F33=DN[:,0,2]*phime[:,0,2]+DN[:,1,2]*phime[:,1,2]+DN[:,2,2]*phime[:,2,2]+DN[:,3,2]*phime[:,3,2]
    
    
    F=np.zeros((me.shape[0],3,3))
    
    F[:,0,0]=F11[:]
    F[:,0,1]=F12[:]
    F[:,0,2]=F13[:]
    
    F[:,1,0]=F21[:]
    F[:,1,1]=F22[:]
    F[:,1,2]=F23[:]
    
    F[:,2,0]=F31[:]
    F[:,2,1]=F32[:]
    F[:,2,2]=F33[:]
    
    Ft=np.zeros((me.shape[0],3,3))
    
    Ft[:,0,0]=F11[:]
    Ft[:,0,1]=F21[:]
    Ft[:,0,2]=F31[:]
    
    Ft[:,1,0]=F12[:]
    Ft[:,1,1]=F22[:]
    Ft[:,1,2]=F32[:]
    
    Ft[:,2,0]=F13[:]
    Ft[:,2,1]=F23[:]
    Ft[:,2,2]=F33[:]
    
    
    Ig=me.reshape((me.shape[0]*4,))
    Mg=np.array([ve0[:],ve0[:],ve0[:],ve0[:]])*1./4.
    Mg=Mg.T
    Mg=Mg.reshape((me.shape[0]*4,))
    M=coo_matrix((Mg,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    Mlumped=np.array(M.diagonal())
    
    
    e11=np.array([ve0[:]*F11[:],ve0[:]*F11[:],ve0[:]*F11[:],ve0[:]*F11[:]])*1./4.
    e11=e11.T
    e11=e11.reshape((me.shape[0]*4,))
    E11=coo_matrix((e11,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E11=np.array(E11.diagonal())
    
    e12=np.array([ve0[:]*F12[:],ve0[:]*F12[:],ve0[:]*F12[:],ve0[:]*F12[:]])*1./4.
    e12=e12.T
    e12=e12.reshape((me.shape[0]*4,))
    E12=coo_matrix((e12,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E12=np.array(E12.diagonal())
    
    
    e13=np.array([ve0[:]*F13[:],ve0[:]*F13[:],ve0[:]*F13[:],ve0[:]*F13[:]])*1./4.
    e13=e13.T
    e13=e13.reshape((me.shape[0]*4,))
    E13=coo_matrix((e13,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E13=np.array(E13.diagonal())
    
    e21=np.array([ve0[:]*F21[:],ve0[:]*F21[:],ve0[:]*F21[:],ve0[:]*F21[:]])*1./4.
    e21=e21.T
    e21=e21.reshape((me.shape[0]*4,))
    E21=coo_matrix((e21,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E21=np.array(E21.diagonal())
    
    
    e22=np.array([ve0[:]*F22[:],ve0[:]*F22[:],ve0[:]*F22[:],ve0[:]*F22[:]])*1./4.
    e22=e22.T
    e22=e22.reshape((me.shape[0]*4,))
    E22=coo_matrix((e22,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E22=np.array(E22.diagonal())
    
    e23=np.array([ve0[:]*F23[:],ve0[:]*F23[:],ve0[:]*F23[:],ve0[:]*F23[:]])*1./4.
    e23=e23.T
    e23=e23.reshape((me.shape[0]*4,))
    E23=coo_matrix((e23,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E23=np.array(E23.diagonal())
    
    e31=np.array([ve0[:]*F31[:],ve0[:]*F31[:],ve0[:]*F31[:],ve0[:]*F31[:]])*1./4.
    e31=e31.T
    e31=e31.reshape((me.shape[0]*4,))
    E31=coo_matrix((e31,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E31=np.array(E31.diagonal())
    
    e32=np.array([ve0[:]*F32[:],ve0[:]*F32[:],ve0[:]*F32[:],ve0[:]*F32[:]])*1./4.
    e32=e32.T
    e32=e32.reshape((me.shape[0]*4,))
    E32=coo_matrix((e32,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E32=np.array(E32.diagonal())
    
    e33=np.array([ve0[:]*F33[:],ve0[:]*F33[:],ve0[:]*F33[:],ve0[:]*F33[:]])*1./4.
    e33=e33.T
    e33=e33.reshape((me.shape[0]*4,))
    E33=coo_matrix((e33,(Ig,Ig)), shape=(q.shape[0],q.shape[0]))
    E33=np.array(E33.diagonal())
    
    P=np.array([E11,E12,E13,E21,E22,E23,E31,E32,E33]).T
    
    return Mlumped, P, ve, ve0


def Write_VTK(xyz, IEN, filename_VTK, scalar_fields=[], scalar_labels=[], vector_fields=[], vector_labels=[]):
    """
    Exports tetrahedral mesh and its data into a VTK file (.vtu)
    xyz           : np.array with nodal coordinates
    IEN           : element connectivity np.array (nodes indexed from 0)
    file_name_VTK : string with file name (including .vtk extension)
    scalar_fields : list with scalar fields in format [data] where,
                        data = np.array with nodal values
    scalar_labels :list with scalar names in format [name] where,
                        name = string with field name
    vector_fields : list with vector fields in format [data] where,
                        data = np.array with vector values as rows
    vector_labels : list with vector names in format [name] where,
                        name = string with vector field name
    Output:
    no output generated, other than the VTK file
    """
    #### Write VTK ####

    cell = np.zeros([IEN.shape[0],IEN.shape[1]+1])
    cell[:,0]  = cell[:,0]+4.0
    cell[:,1:] = IEN
    cell_types = np.zeros([IEN.shape[0],1])+10.0
    f = open(filename_VTK,"w") 
    f.write("# vtk DataFile Version 2.0\n")
    f.write("Unstructured Grid example\n")
    f.write("ASCII\n")
    f.write("DATASET UNSTRUCTURED_GRID\n")
    f.write("POINTS %d float\n" % (xyz.shape[0]))
    np.savetxt(f,xyz,fmt='%-7.6f')
    f.write("CELLS %d %d\n" % (IEN.shape[0],IEN.shape[0]*IEN.shape[1]+IEN.shape[0]))
    np.savetxt(f,cell,fmt='%d')
    f.write("CELL_TYPES %d\n" % (IEN.shape[0]))
    np.savetxt(f,cell_types,fmt='%d')
    f.write("POINT_DATA %d\n" % (xyz.shape[0]))

    # VECTORS:
    for v in np.arange(len(vector_fields)):
        f.write("VECTORS "+vector_labels[v]+" float\n")
        np.savetxt(f,vector_fields[v],fmt='%f')
    
    # SCALARS:
    for s in np.arange(len(scalar_fields)): 
        f.write("SCALARS "+scalar_labels[s]+" float 1\n")
        f.write("LOOKUP_TABLE default\n")
        np.savetxt(f,scalar_fields[s],fmt='%f')
    


def ComputeNodalDefGradTensor(Mlumped,P):
    
    """
    This subroutine computes the nodal values of deformation gradient
    by solving the Mlumped*u=P system, where
    Mlumped    : Diagonal mass matrix (in vector form)
    P          : weighted deformation-gradient vector

    Output:
    Fnodal     : [..., F_nodoi, ...]
    """
    Fnodal = []
    append=Fnodal.append
    nodenumber = 0
    for m, fvec in zip(Mlumped, P[:]) :
        Fnode = np.matrix(np.zeros((3,3)))
        assert (m > 0), "nodal mass is zero in node " + str(nodenumber)
        fvec/=m
        for i in np.arange(3) :
            for J in np.arange(3) :
                Fnode[i,J] = fvec[3*i+J]
        append(Fnode)
        nodenumber += 1
    return Fnodal
    



def Biomech_Analysis(filename_npz,filename_VTK):
    '''
    Estimate deformation measurements and save information in a VTK file and 
    overwrites the info in the .npz file.
    
    Inputs:
        filename_npz  : .npz file containing nodal information
    Outputs:
        file_name_VTK : string with file name (including .vtk extension)
    '''
    
    #print("Loading npz")
    npzfile = np.load(filename_npz)
    LM  = npzfile['elem']
    xyz = npzfile['xyz']
    phi = npzfile['phi']
    
    i_ee = npzfile['intensity_ee']
    i_ei = npzfile['intensity_ei']
    dGF = (i_ei - i_ee)/(-1000.)
    #print("Tet4DefGradSmoothingOptiSPARSE")
    Mlumped, P, ve, ve0 = Tet4DefGradSmoothingOptiSPARSE(xyz, phi, LM)
    #print("ComputeNodalDefGradTensor")
    Fnodal = ComputeNodalDefGradTensor(Mlumped,P)
    M = np.copy(Mlumped)
    Fnodal_T = np.transpose(Fnodal,(0,2,1))
    C        = matrix_multiply(Fnodal_T,Fnodal)

    lam2, N = np.linalg.eig(C[:])    
    for n in np.arange(lam2.shape[0]):
        idx       = np.argsort(lam2[n])[::-1]
        lam2[n,:] = lam2[n,idx]
        N[n]      = N[n,:,idx].T
        
    B = matrix_multiply(Fnodal,Fnodal_T)
    
    aux, nvec = np.linalg.eig(B[:])
    for n in np.arange(aux.shape[0]):
        idx      = np.argsort(aux[n])[::-1]
        aux[n,:] = aux[n,idx]
        nvec[n]  = nvec[n,:,idx].T
    
    lam = lam2**0.5
    J   = np.linalg.det(Fnodal[:])
    #print("Invariants")
    #### Invariantes tensor C ####   
    I1_C = np.sum(lam2,axis=(1))
    I2_C = lam2[:,0]*lam2[:,1]+lam2[:,1]*lam2[:,2]+lam2[:,2]*lam2[:,0]
    I3_C = J**2
    
    #### Invariantes tensor U (I1 e I2 normalizadas)####    
    I1_U = np.sum(lam,axis=(1)) /3.
    I2_U = (lam[:,0]*lam[:,1]+lam[:,1]*lam[:,2]+lam[:,2]*lam[:,0] ) /3.
    I3_U = lam[:,0]*lam[:,1]*lam[:,2]           #### => Jacobian
    
    
    I1_U = np.abs(I1_U)            
    I2_U = np.abs(I2_U)
    I3_U = np.abs(I3_U)
    
    VS = (J-1)*100
    U  = phi - xyz
    #print("ADI, SRI, Others")
    ADI =  ( ((lam[:,0]-lam[:,1])/(lam[:,1]))**2 + ((lam[:,1]-lam[:,2])/(lam[:,2]))**2 )**0.5   

    SRI = ( np.arctan( (lam[:,2]*(lam[:,0]-lam[:,1])) / (lam[:,1]*(lam[:,1]-lam[:,2])) ) ) / (np.pi/2)

    I1_mean = np.sum(M*np.log(I1_U))/np.sum(M) 
    I1_std  = (np.sum(M*(np.log(I1_U)-I1_mean)**2)/np.sum(M))**0.5
    I1_U_LogNorm = (np.log(I1_U) - I1_mean) / I1_std                        
     
    I2_mean = np.sum(M*np.log(I2_U))/np.sum(M) 
    I2_std  = (np.sum(M*(np.log(I2_U)-I2_mean)**2)/np.sum(M))**0.5
    I2_U_LogNorm = (np.log(I2_U) - I2_mean) / I2_std                        

    I3_mean = np.sum(M*np.log(I3_U))/np.sum(M) 
    I3_std  = (np.sum(M*(np.log(I3_U)-I3_mean)**2)/np.sum(M))**0.5
    I3_U_LogNorm = (np.log(I3_U) - I3_mean) / I3_std                        
    
    #print("Saving npz")
    
    np.savez(filename_npz,
             xyz = xyz, elem = LM, phi = phi, DefMap = U, Fnodal = Fnodal, 
             ve0 = ve0, ve = ve, 
             I1_C = I1_C, I2_C = I2_C, I3_C = I3_C, N = N, n = nvec, 
             M = Mlumped, VS = VS,
             I1_U = I1_U, I2_U = I2_U, I3_U = I3_U,
             ADI = ADI, SRI = SRI,
             I1_U_LogNorm = I1_U_LogNorm, 
             I2_U_LogNorm = I2_U_LogNorm, 
             I3_U_LogNorm = I3_U_LogNorm,
             Intensity_Exp = i_ee, 
             Intensity_Insp = i_ei, 
             Delta_Gas_Fraction = dGF)
    
    #print("Defining big lists")
    
    Var_Vectores         = [xyz, phi, U, N[:,:,0], N[:,:,1], N[:,:,2]]
    Var_Vectores_Legend  = ['xyz', 'DefMap', 'DispField', 'N1', 'N2', 'N3']
        
    Var_Escalares        = [ I1_C, I2_C, I3_C,VS, I1_U, I2_U, I3_U, 
                            I1_U_LogNorm, I2_U_LogNorm, I3_U_LogNorm, ADI, SRI,
                            i_ee, i_ei, dGF]
    
    Var_Escalares_Legend = ['Inv_1_C', 'Inv_2_C', 'Inv_3_C','VolStrain', 
                            'Inv_1_U', 'Inv_2_U', 'Inv_3_U', 
                            'Inv_1_U_LogNorm', 'Inv_2_U_LogNorm', 'Inv_3_U_LogNorm', 
                            'ADI', 'SRI', 
                            "Intensity_Exp", "Intensity_Insp", "d_GasFraction"]
    #print("Writing VTK")
    #print("   >  VTK path: %s"%filename_VTK)
    
    #if os.path.isfile(filename_VTK):
    #    print("VTK found before 'Write_VTK' function")
    #else:
    #    print("VTK not found before 'Write_VTK' function")

    Write_VTK(xyz, LM, filename_VTK, Var_Escalares, Var_Escalares_Legend, 
              Var_Vectores, Var_Vectores_Legend)
    
    #if os.path.isfile(filename_VTK):
    #    print("1")
    #else:
    #    print("0")


# =============================================================================
# CODE EXECUTION
# =============================================================================

option = sys.argv[1]
folder = sys.argv[2]

print("*"*80)
print("Program: %s"%sys.argv[0])
print("*"*80+"\n")

print(" * Folder: %s"%folder)
print(" * Option: %s"%option)
print()

if option not in ["1","2","3","4"]:
    print("Option '%s' has no associated program."%option)
    print("WARNING: SHUTTING DOWN PROGRAM")
    exit()

if not os.path.isdir(folder):
    print("Folder '%s' not found"%folder)
    print("WARNING: SHUTTING DOWN PROGRAM")
    exit()
else:
    print("Folder '%s' successfully found!"%folder)

print()

# =============================================================================
# TRANSFORM NII TO INR
# =============================================================================

if option == "1":
    
    print("\nTransforming .nii.gz or .nii to an .inr file")
    
    print(" > Target folder: %s"%folder)
    
    segmentation_image = "NIFTI/NEW_Mask_Exp.nii.gz"
    inr_name = "MESH/geometry.inr"


    if not os.path.isfile(os.path.join(folder,segmentation_image)):
        print("** Segmentation image not found **")
        print("FATAL ERROR")
        exit()
    
    # =========================================================================
    # Reading the segmentation image
    # =========================================================================
    # Open the image in order to extract relevant dimension data that will also
    # be useful for later computations.
    print("Inspecting segmentation")
    seg = nib.load(os.path.join(folder,segmentation_image))
    print("   + Extracting data")
    seg_data = np.array(seg.dataobj)
    sx, sy, sz = seg_data.shape
    dx, dy, dz, _ = np.diagonal(seg.affine).__abs__()
    print("   + Done")
    
    # =========================================================================
    # Writing an associated INR file
    # =========================================================================
    # INR is a poorly documentated filetype associated with CGAL. Documentation
    # regarding it's structure is oddly found and it's mostly secondary sources
    # (people commenting on threads). Currently it's working fine after being
    # translated from Python 2.7.
    
    # write INR
    print(" > Writing INR")
    writeINR(os.path.join(folder,inr_name), seg_data, dx, dy, dz)
    print("   + Done")
    
# =============================================================================
# TRANSFORM MESH INTO AN .NPZ FILE
# =============================================================================
elif option == "2":

    ''' processing the .mesh CGAL output file '''
    
    print("Transforming the .mesh into a ready-to-be-used FE mesh")
    
    # mesh name
    name = "triang.mesh"
    # read the original mesh file
    
    try:
        print("\n > Reading '%s' file:"%name.split("/")[-1])
        mesh = io.read(os.path.join(folder,name))
        print("   + Success! ")
    except:
        print("   + Failure")
        exit()
        
    faces = mesh.cells_dict['triangle']
    
    # generate a trimesh object
    trim = trimesh.Trimesh(vertices=mesh.points,
                              faces=mesh.cells_dict['triangle'])
    
    try:
        # smooth
        print("\n > Smoothing")
        trimesh.smoothing.filter_humphrey(trim) 
        print("   + Success!")
        # save the smoothed mesh to vtk
        io.write_points_cells(os.path.join(folder,"surface.vtk"),
                              trim.vertices,
                              {"triangle":trim.faces})
        print("   + Saving as surface.vtk")
        
    except:
        print("   + Failure")
        exit()
        

    
    # tetrahedralize surface mesh using tetgen (couldn't run the cgal routine)
    # tetgen has too many parameters. this could be improved.
    
    try:
        print(" > Generating tetrahedral mesh")
        trim_lung = tetgen.TetGen(trim.vertices, trim.faces)
        tet_lung = trim_lung.tetrahedralize(order=1, 
                                        minratio = 1.1,
                                        quality=True)
        print("   + Success!")
    except:
        print("   + Failure")
        exit()

    xyz = tet_lung[0]
    IEN = tet_lung[1]
    
    try:
        print(" > Exaluating flying nodes")
        xyz, IEN = Evaluate_FlyingNodes(xyz, IEN)
        print("   + Success!")
    except:
        print("   + Failure")
        exit()
    
    np.savez(os.path.join(folder,"surface_new.npz"), xyz=xyz, elem=IEN)
    #io.write_points_cells(os.path.join(folder, "MESH/Exp.vtk"), 
    #                      points=xyz, cells={"tetra":IEN})
    
    #if os.path.isfile(os.path.join(folder, "MESH/Exp.vtk")):
    #    print("\n > Exp.vtk successfully generated")

    if os.path.isfile(os.path.join(folder, "surface_new.npz")):
        print(" > Exp_NEW.npz successfully generated")
        
    if os.path.isfile(os.path.join(folder, "MESH/geometry.inr")):
        os.remove(os.path.join(folder, "MESH/geometry.inr"))
        print(" > Removed the geometry.inr file successfully")
        
if option == "3":
    
    ''' interpotation of intensities and the displacement field '''
    
    print("\nInterpolating displacements and intensities")
    # 1) check for the cpp file and the target images
    if not os.path.isfile(os.path.join(folder, 
                                   "REGISTRATION/cpp_9000-000666.nii.gz")):
        print(" > CPP file not found!")
        print(" > SHUTTING DOWN PROGRAM")
        exit()
    else:
        print(" > CPP file found")
    
    for img in ["Exp", "Insp"]:
        if not os.path.isfile(os.path.join(folder, 
                                       "NIFTI/%s.nii.gz"%img)):
            print(" > %s.nii.gz file not found!"%img)
            print(" > SHUTTING DOWN PROGRAM")
            exit()
        else:
            print(" > %s.nii.gz file found"%img)
        
    # 2) check for the prepared mesh
    if not os.path.isfile(os.path.join(folder, "MESH/Exp_NEW.npz")):
        print(" > NEW_Exp.npz not found")
        print(" > SHUTTING DOWN PROGRAM")
        exit()
    else:
        print(" > NEW_Exp.npz found")

    
    # 3) interpolate the displacement field
    print("\nInterpolating displacements")
    try:
        InterpolateBSplines(os.path.join(folder, "MESH/Exp_NEW.npz"), 
                        os.path.join(folder, "REGISTRATION/cpp_9000-000666.nii.gz"), 
                        os.path.join(folder, "NIFTI/Exp.nii.gz"))
        print(" > Success")
    except:
        print(" > Failure")
        exit()

    # 4) interpolate the intensities from the images
    npz = os.path.join(folder, "MESH/Exp_NEW.npz")
    img_ref = os.path.join(folder, "NIFTI/Exp.nii.gz")
    img_flo = os.path.join(folder, "NIFTI/Insp.nii.gz")
    
    print("\nInterpolating intensities")
    try:
        i_ee, i_ei = interpolateIntensities(npz, img_ref, img_flo)
        print(" > Success")
    except:
        print(" > Failure")
        exit()

    
    # 5) save the data
    print("\nSaving results")

    # load the mesh
    msh = np.load(npz)
    # update the mesh
    np.savez(npz, xyz=msh['xyz'], elem=msh['elem'], phi=msh['phi'],
             intensity_ee=i_ee, intensity_ei=i_ei)
    
    # save as a vtk object
    #io.write_points_cells(os.path.join(folder, "MESH/Exp.vtk"), 
    #                      points=msh['xyz'], 
    #                      point_data={'defmap':msh['phi'],
    #                                  'dispfield':(msh['phi']-msh['xyz']),
    #                                  'intensity_ee':i_ee,
    #                                  'intensity_ei':i_ei,
    #                                  'delta_gas_fraction':(i_ei-i_ee)/-1000}, 
    #                      cells={"tetra":msh['elem']})
    print(" > Done")
    
# =============================================================================
# Biomechanical analysis
# =============================================================================

elif option == "4":
    
    print("Option 4: Biomechanical analysis \n")
   

    npz = os.path.join(folder, "MESH/Exp_NEW.npz")

    if not os.path.isfile(os.path.join(folder, "MESH/Exp_NEW.npz")):
        print(" > NEW_Exp.npz not found")
        print(" > SHUTTING DOWN PROGRAM")
        exit()
    else:
        print(" > NEW_Exp.npz found")
        
        
    vtk = os.path.join(folder, "MESH/Exp.vtk")

    if os.path.isfile(vtk):
        os.remove(vtk)

    try:  
        Biomech_Analysis(npz, vtk)
        print(" > Success!")
    except: 
        print(" > Failure")
        exit()
