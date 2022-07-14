'''Model file for alveolar perfusion and gas exchange simulations used in
Zurita & Hurtado, 2022.
'''

__author__ = 'pzuritas'
__email__ = 'pzurita@uc.cl'

import os
import time as tm
from datetime import date
import numpy as np
from fenics import *
from dolfin import *
from src.boundaries import *
from ufl.algorithms.formtransformations import compute_form_arities


class PerfusionGasExchangeModel():
    '''FEniCS simulater class for microscale alveolar perfusion and gas exchange
    model.
    '''
    def __init__(self, folder_path, params):
        '''Instance the model.

        folder_path: path to folder for storing solution files. (string)
        params: dictionary for parameter values. (dict)
        '''
        self.folder_path = folder_path
        self.params = params

        # Save details
        if not os.path.exists(self.folder_path):
            os.mkdir(self.folder_path)
        with open(self.folder_path+'/info.txt', 'w') as file:
            file.write(f'Simulation done on {date.today()}\n\n')
            file.write(f'Parameters used:\n\n')
            for param in self.params:
                file.write(f'Parameter {param}: {self.params[param]}\n')

    def import_mesh(self, mesh_path, meshtype=None, type="h5", periodic=False):
        '''Imports mesh from .h5 file for use in simulations.

        mesh_path: path to file. (string)
        meshtype: specific case for refinement. (string)
        type: file format for mesh. (string)
        periodic: use periodic boundary conditions. (bool)
        '''
        if type == "h5":
            hdf5 = HDF5File(MPI.comm_world, mesh_path, 'r')
            self.mesh = Mesh()
            hdf5.read(self.mesh, 'mesh', False)
        elif type == "xml":
            self.mesh = Mesh(mesh_path)
        elif type == "xdmf":
            self.mesh = Mesh()
            with XDMFFile(mesh_path) as infile:
                infile.read(self.mesh)
        else:
            raise ValueError("type of mesh should be h5 or xml")
        
        print("Mesh imported")
        dir_arr_flow = np.array(
                [coords[0] for coords in self.mesh.coordinates()]
        )
        dir_arr_rest = np.array(
                [coords[1] for coords in self.mesh.coordinates()]
        )
        # Node coordinates for principal (x) direction
        self.dir_max_flow = np.max(dir_arr_flow)  # Maximum principal direction coordinate
        self.dir_min_flow = np.min(dir_arr_flow)  # Minimum principal direction coordinate
        self.dir_max = np.max(dir_arr_rest)
        self.dir_min = np.min(dir_arr_rest)
        print("self.dir_max_flow = ", self.dir_max_flow)
        print("self.dir_min_flow = ", self.dir_min_flow)
        print("self.dir_max = ", self.dir_max)
        print("self.dir_min = ", self.dir_min)
        len_flow = self.dir_max_flow - self.dir_min_flow
        len_air = self.dir_max - self.dir_min
        delta = (len_flow - len_air)/2

        # Save mesh dims for other uses

        self.dims = (len_flow, self.dir_max, self.dir_max)

        # Flag for periodicity

        self.periodic = periodic


    def instance_boundaries(self, mesh=None):
        '''Instances the relevant boundaries for boundary conditions.
        
        mesh: type of mesh created. None, "slab" or "tkd". (None or str)
        '''
        
        self.meshtype = mesh
        # Instance the relevant boundaries
        print("Instancing boundaries")
        
        tol = 0.5

        if not self.periodic:
            self.gamma_in = GammaIn(
                self.dir_min_flow, self.dir_max_flow, tol
            )
            self.gamma_out = GammaOut(
                self.dir_min_flow, self.dir_max_flow, tol
            )
            self.gamma_air = GammaAir(
                self.dir_min, self.dir_max, tol
            )

            print("gamma_in.dir_min = ", self.gamma_in.dir_min)
            print("gamma_in.dir_max = ", self.gamma_in.dir_max)
            print("gamma_out.dir_min = ", self.gamma_out.dir_min)
            print("gamma_out.dir_max = ", self.gamma_out.dir_max)

            print("gamma_air.dir_min = ", self.gamma_air.dir_min)
            print("gamma_air.dir_max = ", self.gamma_air.dir_max)

            # Declare the boundaries in the mesh and tag them

            self.boundaries = MeshFunction('size_t', self.mesh, dim=2)
            self.boundaries.set_all(3)
            self.gamma_in.mark(self.boundaries, 1)
            self.gamma_out.mark(self.boundaries, 2)
            self.gamma_air.mark(self.boundaries, 3)

        boundaries = File(self.folder_path+'/bnd/bnd.pvd')
        boundaries << self.boundaries

    def instance_function_spaces(self):
        '''Instances the relevant function spaces.'''

        if not self.periodic:
            self.W_h = FunctionSpace(self.mesh, 'Lagrange', 2) # degree n+1
            self.V_h = VectorFunctionSpace(self.mesh, 'Lagrange', 1) # degree n
        else:
            self.W_h = FunctionSpace(
                self.mesh, 'Lagrange', 2,
                constrained_domain=self.gamma_pi
            )
            self.V_h = VectorFunctionSpace(
                self.mesh, 'Lagrange', 1,
                constrained_domain=self.gamma_pi
            )

    def sim_p(self, save=True, meshtype=None):
        '''Solves the perfusion (P) problem of the model.
        
        save: saves to vtk. (bool)
        meshtype: type of mesh. None, "slab" or "tkd". (None or str)
        '''
        print("began sim_p")
        self.instance_boundaries(mesh=meshtype)
        print("boundaries instanced")
        self.instance_function_spaces()

        # Declare Dirichlet boundary conditions for (P)

        self.p_dbc = [
            DirichletBC(self.W_h, self.params['p_min'], self.gamma_out)
        ]
        # Assemble problem

        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        p = TrialFunction(self.W_h)
        v = TestFunction(self.W_h)
        a = inner(grad(p), grad(v))*dx
        F = Constant(0)*v*dx
        F += self.params["u_in"]*self.params["mu"]/self.params["kappa"]*v*ds(1)

        # Solve problem

        self.p = Function(self.W_h)
        solve(
            a == F, self.p, self.p_dbc, solver_parameters={'linear_solver':'mumps'}#,
            #solver_parameters={
            #    'linear_solver': 'gmres',
            #    'preconditioner': 'ilu'
            #}
        )
        print("P problem solved")
        self.u = project(
            -1/self.params['mu']*self.params['kappa']*grad(self.p),
            self.V_h
        )
        print("u problem solved")
        self.u.rename("u", "blood velocity [um/s]")
        self.p.rename("p", "blood pressure [mmHg]")

        if save:

            # Save solutions
            print("Saving solutions")
            u_file = File(self.folder_path+'/p/u_all_flat_1_fixed.pvd')
            u_file << self.u
            p_file = File(self.folder_path+'/p/p_all_flat_1_fixed.pvd')
            p_file << self.p
            print("saved")

    def set_u(self, value=(0, 0, 0), meshtype='slab', save=True):
        '''Prescribes a velocity field u to the mesh instead of solving (P).

        value: uniform velocity field. (tuple)
        meshtype: specific case of mesh. (string or None)
        save: saves to vtk. (bool)
        '''
        self.instance_boundaries(mesh=meshtype)
        self.instance_function_spaces()

        self.u = project(Expression(tuple(map(str, value)), degree=1), self.V_h)
        self.u.rename("u", "Blood velocity")

        if save:

            u_file = File(self.folder_path+'/p/u.pvd')
            u_file << self.u

    def f(self, X, p_X, c_HbX, c_HbY):
        '''Generation rate of X as defined in the model.
        
        X: gas species. (string)
        p_X: partial pressure of X. (FEniCS Function)
        c_HbX: concentration of Hb(X). (FEniCS Function)
        c_HbY: concentration of Hb(Y). (FEniCS Function)
        '''
        c_t = self.params['c_t']
        if X == 'O2':
            beta_O2 = self.params['beta_O2']
            k_O2 = self.params['k_O2']
            k_prime_O2 = self.params['k_prime_O2']
            first_term = k_prime_O2*(c_t - c_HbX - c_HbY)*p_X
            alpha = 1
            second_term = -k_O2/beta_O2*c_HbX
            beta = 1
            return alpha*first_term + beta*second_term
        elif X == 'CO2':
            beta_CO2 = self.params['beta_CO2']
            k_CO2 = self.params['k_CO2']
            k_prime_CO2 = self.params['k_prime_CO2']
            first_term = k_prime_CO2*(c_t - c_HbX - c_HbY)*p_X
            alpha = 1
            second_term = -k_CO2/beta_CO2*c_HbX
            beta = 1
            return alpha*first_term + beta*second_term
        else:
            raise ValueError('Gas species in f must be O2 or CO2.')

    def g(self, X, p_X, c_HbX, c_HbY):
        '''Scaled generation rate of X as defined in the model.
        
        X: gas species. (string)
        p_X: partial pressure of X. (FEniCS Function)
        c_HbX: concentration of Hb(X). (FEniCS Function)
        c_HbY: concentration of Hb(Y). (FEniCS Function)
        '''
        if X == 'O2':
            beta_O2 = self.params['beta_O2']
            return beta_O2*self.f(X, p_X, c_HbX, c_HbY)
        elif X == 'CO2':
            beta_CO2 = self.params['beta_CO2']
            return beta_CO2*self.f(X, p_X, c_HbX, c_HbY)
        else:
            raise ValueError('Gas species in f must be O2 or CO2.')

    def sim_sbst(self, hb=True, save=True, guess=None, stabilizer=None):
        '''Solves the steady state blood-side transport (T) problem of the
        model.
        
        hb: toggle effects of hemoglobin, that is, chooses between solving
        (WT) or (WLT). (bool)
        save: saves to vtk. (bool)
        guess: starting point for Newton iterations. s^0_L in the paper.
        (FEniCS Function or None)
        '''
        # Instance parameters

        p_air_O2 = self.params['p_air_O2']
        d_ba_O2 = self.params['d_ba_O2']
        d_pla_O2 = self.params['d_pla_O2']

        p_air_CO2 = self.params['p_air_CO2']
        d_ba_CO2 = self.params['d_ba_CO2']
        d_pla_CO2 = self.params['d_pla_CO2']

        h_ba = self.params['h_ba']

        # Instance function space for the multi-field problem

        element = VectorElement('P', tetrahedron, 1, dim=4)
        self.M_h = FunctionSpace(self.mesh, element)
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        # Declare functions and test functions, and initial point for iterations

        if not guess:
            x = Function(self.M_h)
        else:
            x = project(guess, self.M_h)
        self.p_O2, self.p_CO2, self.c_HbO2, self.c_HbCO2 = split(x)

        v, w, eta, xi = TestFunctions(self.M_h)

        n = FacetNormal(self.mesh)

        # Define residuals

        G_p_O2 = d_pla_O2*inner(grad(self.p_O2), grad(v))*dx
        G_p_O2 += -d_ba_O2/h_ba*Constant(p_air_O2)*v*ds(3)
        G_p_O2 += d_ba_O2/h_ba*self.p_O2*v*ds(3)
        G_p_O2 += inner(self.p_O2*self.u, n)*v*ds(2)
        G_p_O2 += -inner(self.p_O2*self.u, grad(v))*dx
        if hb:
            G_p_O2 += self.f('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*v*dx

        G_p_CO2 = d_pla_CO2*inner(grad(self.p_CO2), grad(w))*dx
        G_p_CO2 += -d_ba_CO2/h_ba*Constant(p_air_CO2)*w*ds(3)
        G_p_CO2 += d_ba_CO2/h_ba*self.p_CO2*w*ds(3)
        G_p_CO2 += inner(self.p_CO2*self.u, n)*w*ds(2)
        G_p_CO2 += -inner(self.p_CO2*self.u, grad(w))*dx
        if hb:
            G_p_CO2 += self.f('CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2)*w*dx

        if hb:
            G_c_O2 = self.g('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*eta*dx
            G_c_O2 += inner(self.c_HbO2*self.u, grad(eta))*dx
            G_c_O2 += -inner(self.c_HbO2*self.u, n)*eta*ds(2)

            G_c_CO2 = self.g(
                'CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2
            )*xi*dx
            G_c_CO2 += inner(self.c_HbCO2*self.u, grad(xi))*dx
            G_c_CO2 += -inner(self.c_HbCO2*self.u, n)*xi*ds(2)
        else:
            G_c_O2 = self.c_HbO2*eta*dx
            G_c_CO2 = self.c_HbCO2*xi*dx
        
        if guess:
            if stabilizer == 'SUPG':
               
                # Norm of velocity field and element size
                norm_u = norm(self.u)
                h = self.mesh.hmax()
                
                # Stabilizing parameters (tau), with approximation
                
                # p_O2 stabilizing parameter
                k = d_pla_O2
                Pe = (norm_u*h)/(2*k) # value =~ 24
                if Pe >= 1:
                    tau_p_O2 = h/(2*norm_u)
                else: 
                    tau_p_O2 = (h**2)/(12*k)
                
                # p_CO2 stabilizing parameter
                k = d_pla_CO2
                Pe = (norm_u*h)/(2*k)
                if Pe >= 1:
                    tau_p_CO2 = h/(2*norm_u)
                else: 
                    tau_p_CO2 = (h**2)/(12*k)
                    
                # Temporary c_HbO2 and c_HbCO2 stabilizing parameters 
                tau_c_HbO2 = tau_p_O2*1e-16
                tau_c_HbCO2 = tau_p_CO2*1e-16

                # Advection terms 
                L_adv_p_O2 =    inner(self.u,grad(v))
                L_adv_p_CO2 =   inner(self.u, grad(w))
                L_adv_c_HbO2 =  -inner(self.u, grad(eta))
                L_adv_c_HbCO2 = -inner(self.u, grad(xi))
                
                # LHS of (47a) and (47b) equations
                L_p_O2 =        -d_pla_O2*div(grad(self.p_O2))      + inner(self.u,grad(self.p_O2))     + self.f('O2', self.p_O2, self.c_HbO2, self.c_HbCO2) #
                L_p_CO2 =       -d_pla_CO2*div(grad(self.p_CO2))    + inner(self.u,grad(self.p_CO2))    + self.f('CO2', self.p_CO2, self.c_HbO2, self.c_HbCO2) #
                L_c_HbO2 =      self.g('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)   - inner(self.u,grad(self.c_HbO2)) 
                L_c_HbCO2 =     self.g('CO2', self.p_CO2, self.c_HbO2, self.c_HbCO2) - inner(self.u,grad(self.c_HbCO2)) 

                G_p_O2  += tau_p_O2*L_p_O2*L_adv_p_O2*dx
                G_p_CO2 += tau_p_CO2*L_p_CO2*L_adv_p_CO2*dx
                G_c_O2  += tau_c_HbO2*L_adv_c_HbO2*L_c_HbO2*dx
                G_c_CO2 += tau_c_HbCO2*L_adv_c_HbCO2*L_c_HbCO2*dx
                print("Added all SUPG stabilizing terms.")

            elif stabilizer == 'GLS':
                tau_k = 1
                # Norm of velocity field and element size
                norm_u = norm(self.u)
                h = self.mesh.hmax()
                
                # Stabilizing parameters (tau), with approximation
                
                # p_O2 stabilizing parameter
                k = d_pla_O2
                
                L_v =           -d_pla_O2*div(grad(v)) + inner(self.u,grad(v)) + self.f('O2',  v, self.c_HbO2, self.c_HbCO2)
                L_w =           -d_pla_CO2*div(grad(w))+ inner(self.u,grad(w)) - self.f('CO2', w, self.c_HbO2, self.c_HbCO2)
                L_eta =         self.g('O2', eta, self.c_HbO2, self.c_HbCO2) - inner(self.u,grad(eta))
                L_xi =          self.g('CO2', xi, self.c_HbO2, self.c_HbCO2) - inner(self.u,grad(xi))

                L_p_O2 =        -d_pla_O2*div(grad(self.p_O2))  + inner(self.u,grad(self.p_O2))  + self.f('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)
                L_p_CO2 =       -d_pla_CO2*div(grad(self.p_CO2))+ inner(self.u,grad(self.p_CO2)) - self.f('CO2', self.p_CO2, self.c_HbO2, self.c_HbCO2)
                L_c_HbO2 =      self.g('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)   - inner(self.u,grad(self.c_HbO2))
                L_c_HbCO2 =     self.g('CO2', self.p_CO2, self.c_HbO2, self.c_HbCO2) - inner(self.u,grad(self.c_HbCO2)) 

                G_p_O2 += tau_k*inner(L_v,L_p_O2)*dx
                G_p_CO2 += tau_k*inner(L_w,L_p_CO2)*dx
                G_c_O2 += tau_k*inner(L_eta,L_c_HbO2)*dx
                G_c_CO2 += tau_k*inner(L_xi,L_c_HbCO2)*dx

            elif stabilizer == 'SGS':
                tau_k = 1
                L_v_star =      -d_pla_O2*div(grad(v)) + dot(self.u,grad(v)) - self.f('O2',  v, self.c_HbO2, self.c_HbCO2)
                L_w_star =      -d_pla_CO2*div(grad(w))+ dot(self.u,grad(w)) - self.f('CO2', w, self.c_HbO2, self.c_HbCO2)
                L_eta_star =    - self.g('O2', eta, self.c_HbO2, self.c_HbCO2) - dot(self.u,grad(eta))
                L_xi_star =     - self.g('CO2', xi, self.c_HbO2, self.c_HbCO2) - dot(self.u,grad(xi))

                L_p_O2 =        -d_pla_O2*div(grad(self.p_O2))+ dot(self.u,grad(self.p_O2)) + self.f('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*v
                L_p_CO2 =       -d_pla_CO2*div(grad(self.p_CO2))+ dot(self.u,grad(self.p_CO2)) - self.f('CO2', self.p_CO2, self.c_HbO2, self.c_HbCO2)*w
                L_c_O2 =        self.g('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*eta - dot(self.u,grad(self.c_HbO2))
                L_c_CO2 =       self.g('CO2', self.p_CO2, self.c_HbO2, self.c_HbCO2)*xi - dot(self.u,grad(self.c_HbCO2)) 

                G_p_O2 += tau_k*inner(L_v_star,L_p_O2)*dx
                G_p_CO2 += tau_k*inner(L_w_star,L_p_CO2)*dx
                G_c_O2 += tau_k*inner(L_eta_star,L_c_HbO2)*dx
                G_c_CO2 += tau_k*inner(L_xi_star,L_c_HbCO2)*dx

        G = G_p_O2 + G_p_CO2 + G_c_O2 + G_c_CO2
     
        if save:
        
            # Create files for output

            p_O2_file = File(self.folder_path+'/sbst/pO2_all_flat_1_fixed.pvd')
            p_CO2_file = File(self.folder_path+'/sbst/pCO2_all_flat_1_fixed.pvd')
            c_HbO2_file = File(self.folder_path+'/sbst/cHbO2_all_flat_1_fixed.pvd')
            c_HbCO2_file = File(self.folder_path+'/sbst/cHbCO2_all_flat_1_fixed.pvd')

        # Declare Dirichlet boundary conditions for (SBST)

        self.sbst_dbc = [
            DirichletBC(
                self.M_h.sub(0), Constant(self.params['p_O2_in']), self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(1), Constant(self.params['p_CO2_in']),
                self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(2), Constant(self.params['c_HbO2_in']),
                self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(3), Constant(self.params['c_HbCO2_in']),
                self.gamma_in
            )
        ]

        # Solve variational problem

        solve(
            G == 0, x, self.sbst_dbc,
            solver_parameters={"newton_solver": {
                "relative_tolerance": 1E-8,
                "absolute_tolerance": 1E-8#,
                #"linear_solver": "gmres",
                #"preconditioner": "ilu"
            }}
        )

        if save:

            # Save solution to files

            _p_O2, _p_CO2, _c_HbO2, _c_HbCO2 = x.split()
            _p_O2.rename("p_O2", "oxygen partial pressure [mmHg]")
            _p_CO2.rename("p_CO2", "carbon dioxide partial pressure [mmHg]")
            _c_HbO2.rename("c_HbO2", "oxyhemoglobin concentration [mole/um^3]")
            _c_HbCO2.rename(
                "c_HbCO2", "carboaminohemoglobin concentration [mole/um^3]"
            )
            p_O2_file << _p_O2
            p_CO2_file << _p_CO2
            c_HbO2_file << _c_HbO2
            c_HbCO2_file << _c_HbCO2

        return x

    def compute_airflow(self):
        '''Experimental. Tries to compute airflow across the blood-air barrier.
        '''
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        n = FacetNormal(self.mesh)
        o2 = assemble(
            dot(grad(self.p_O2), n)*ds(3) + Constant(0)*dx(
                domain=self.mesh, subdomain_data=self.boundaries
            )
        )
        co2 = assemble(
            dot(grad(self.p_CO2), n)*ds(3) + Constant(0)*dx(
                domain=self.mesh, subdomain_data=self.boundaries
            )
        )
        return o2, co2

    def compute_blood_conservation(self):
        '''Experimental. Tries to compute inflow vs outflow of blood.
        '''
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)
        n = FacetNormal(self.mesh)
        inlet = assemble(
            dot(self.u, n)*ds(1) + Constant(0)*dx(
                domain=self.mesh, subdomain_data=self.boundaries
            )
        )
        outlet = assemble(
            dot(self.u, n)*ds(2) + Constant(0)*dx(
                domain=self.mesh, subdomain_data=self.boundaries
            )
        )
        return inlet + outlet