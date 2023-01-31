'''Boundary classes for alveolar perfusion and gas exchange simulations.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

from dolfin import *
import numpy as np

class GammaIn(SubDomain):
    '''Subdomain class for boundary conditions on the inlet of a rectangular
    prism.
    '''

    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''
        
        return on_boundary and near(x[0], self.dir_min, self.tol)

class GammaOut(SubDomain):
    '''Subdomain class for boundary conditions on the outlet of a rectangular
    prism.
    '''
    def __init__(self, dir_min, dir_max, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min = dir_min
        self.dir_max = dir_max
        self.tol = tol
    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''

        return on_boundary and near(x[0], self.dir_max, self.tol)

class GammaAir(SubDomain):
    '''Subdomain class for boundary conditions on the blood-air barrier of a
    rectangular prism.
    '''
    def __init__(self, dir_min_y, dir_max_y, dir_min_z, dir_max_z, tol):
        '''Instance the subdomain.
        
        dir_min: minimum value of flow direction. (float)
        dir_max: maximum value of flow direction. (float)
        tol: tolerance for numerical roundoff in element tagging. (float)
        '''
        super().__init__()
        self.dir_min_y = dir_min_y
        self.dir_max_y = dir_max_y
        self.dir_min_z = dir_min_z
        self.dir_max_z = dir_max_z        
        self.tol = tol

    def inside(self, x, on_boundary):
        '''Checks if position is on subdomain.
        
        x: position.
        on_boundary: True if on element boundary. (bool)
        '''

        return on_boundary and (near(x[1], self.dir_min_y, self.tol) or near(x[1], self.dir_max_y, self.tol) or near(x[2], self.dir_min_z, self.tol) or near(x[2], self.dir_max_z, self.tol))
    
class InletOutlet(SubDomain):
    def inside(self, x, on_boundary):
        return (x[0]<5-DOLFIN_EPS or x[0]>40-5+DOLFIN_EPS)
    
class Inlet(SubDomain):
    def inside(self, x, on_boundary):
        return x[0]<5-DOLFIN_EPS