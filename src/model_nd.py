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

class NondimensionalPerfusionGasExchangeModel():
    '''FEniCS simulater class for microscale alveolar perfusion and gas exchange
    model.
    '''
    def __init__(self, folder_path, params, solver='', f_dim = 2, vecf_dim = 1):
        '''Instance the model.

        folder_path: path to folder for storing solution files. (string)
        params: dictionary for parameter values. (dict)
        '''
        self.folder_path = folder_path
        self.params = params
        self.solver = solver
        self.f_dim = f_dim
        self.vecf_dim = vecf_dim

        if not os.path.exists(self.folder_path):
            os.makedirs(self.folder_path, exist_ok=True)
        with open(self.folder_path+'/info.txt', 'w') as file:
            file.write(f'Simulation done on {date.today()}\n\n')
            file.write(f'Parameters used:\n\n')
            for param in self.params:
                file.write(f'Parameter {param}: {self.params[param]}\n')

    def instance_nondimensional_params(self):
        U = self.params['u_in']
        L = 6
        rho = 1E-3 # uPa s^2 um^-2
        mu = 3499.7 # uPa s
        
        conv = 1E6*133.322 # uPa / mmHg
        
        self.U = U
        self.L = L
        self.rho = rho
        self.mu = mu
        self.conv = conv
        
        self.Re = (rho*U*L)/mu
        self.Tau = mu*U/L
        self.Pe_O2 = L*U/(self.params['d_pla_O2'])
        self.Pe_CO2 = L*U/(self.params['d_pla_CO2'])
        self.Dap_T_O2 = self.params['k_prime_O2']*self.params['c_t']*(L**2)/self.params['d_pla_O2']
        self.Dap_T_CO2 = self.params['k_prime_CO2']*self.params['c_t']*(L**2)/self.params['d_pla_CO2']
        self.Dap_O2 = self.params['k_prime_O2']*self.params['c_HbO2_in']*(L**2)/self.params['d_pla_O2']
        self.Dap_CO2 = self.params['k_prime_CO2']*self.params['c_HbCO2_in']*(L**2)/self.params['d_pla_CO2']
        self.Da_O2 = self.params['k_O2']*(L**2)/self.params['d_pla_O2']
        self.Da_CO2 = self.params['k_CO2']*(L**2)/self.params['d_pla_CO2']        
        self.H_O2  = conv*self.params['c_HbO2_in']/(rho*(U**2)*self.params['beta_O2'])
        self.H_CO2 = conv*self.params['c_HbCO2_in']/(rho*(U**2)*self.params['beta_CO2'])

        print(f"Re = {self.Re}")
        print(f"Tau = {self.Tau}")
        print(f"Pe_O2 = {self.Pe_O2}")
        print(f"Pe_CO2 = {self.Pe_CO2}")
        print(f"Da'_(T, O2) = {self.Dap_T_O2}")
        print(f"Da'_(T, CO2) = {self.Dap_T_CO2}")
        print(f"Da'_O2 = {self.Dap_O2}")
        print(f"Da'_CO2 = {self.Dap_CO2}")
        print(f"Da_O2 = {self.Da_O2}")
        print(f"Da_CO2 = {self.Da_CO2}")
        print(f"H_O2 = {self.H_O2}")
        print(f"H_CO2 = {self.H_CO2}")

        self.S_O2_T = (mu*self.params['beta_O2']*self.params['k_prime_O2'])*self.params['c_t']/(self.params['c_HbO2_in']) * (1/conv)
        self.S_CO2_T = (mu*self.params['beta_CO2']*self.params['k_prime_CO2'])*self.params['c_t']/(self.params['c_HbCO2_in']) * (1/conv)
        self.S_O2_O2 = mu*self.params['beta_O2']*self.params['k_prime_O2'] * (1/conv)
        self.S_CO2_CO2 = mu*self.params['beta_CO2']*self.params['k_prime_CO2']* (1/conv)
        self.S_O2_CO2 = (mu*self.params['beta_O2']*self.params['k_prime_O2']*self.params['c_HbCO2_in'])/self.params['c_HbO2_in']* (1/conv)
        self.S_CO2_O2 = (mu*self.params['beta_CO2']*self.params['k_prime_CO2']*self.params['c_HbO2_in'])/self.params['c_HbCO2_in']* (1/conv)
        self.W_O2 = L*self.params['k_O2']/U
        self.W_CO2 = L*self.params['k_CO2']/U
        self.R_O2 = (self.params['d_ba_O2']*L)/(self.params['d_pla_O2']*self.params['h_ba'])
        self.R_CO2 = (self.params['d_ba_CO2']*L)/(self.params['d_pla_CO2']*self.params['h_ba'])

        print(f"S_(O2, T) = {self.S_O2_T}")
        print(f"S_(CO2, T) = {self.S_CO2_T}")
        print(f"S_(O2, O2) = {self.S_O2_O2}")
        print(f"S_(CO2, CO2) = {self.S_CO2_CO2}")
        print(f"S_(O2, CO2) = {self.S_O2_CO2}")
        print(f"S_(CO2, O2) = {self.S_CO2_O2}")
        print(f"W_O2 = {self.W_O2}")
        print(f"W_CO2 = {self.W_CO2}")
        print(f"R_O2 = {self.R_O2}")
        print(f"R_CO2 = {self.R_CO2}")

    def import_mesh(self, mesh_path=os.path.join("../raw-data/40_h_repaired", "40_h_repaired.xdmf"), meshtype=None, type="h5", periodic=False, max_dims=[41.2, 40.7, 40.7], min_dims=[-1.2, -0.7, -0.7], tol=0.2, box_side_length=None, box_nodes=None):
        '''Imports mesh from .h5 file for use in simulations.
        mesh_path: path to file. (string)
        meshtype: specific case for refinement. (string)
        type: file format for mesh. (string)
        periodic: use periodic boundary conditions. (bool)
        '''
        
        if box_side_length is None:
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
            
        else:
            if box_nodes is None:
                raise ValueError("Parameter box_nodes_amount should be non-zero.")
            else:
                max_dims = [box_side_length, box_side_length, box_side_length]
                min_dims = [0,0,0]
                self.mesh = BoxMesh(Point(0,0,0), Point(box_side_length,box_side_length,box_side_length), box_nodes, box_nodes, box_nodes)
                print("Created b")

        print("Coordinates have shape ", np.shape(self.mesh.coordinates()))

        self.dir_max_flow =  max_dims[0]
        self.dir_min_flow =  min_dims[0]
        self.dir_max_y =  max_dims[1]
        self.dir_min_y =  min_dims[1]
        self.dir_max_z =  max_dims[2]
        self.dir_min_z =  min_dims[2]

        len_flow = self.dir_max_flow - self.dir_min_flow
        len_y = self.dir_max_y - self.dir_min_y
        len_z = self.dir_max_z - self.dir_min_z
        
        # Tolerance for boundary instancing
        self.tol = tol

        # Save mesh dims for other uses
        self.dims = (len_flow, len_y, len_z)

        # Flag for periodicity
        self.periodic = periodic


    def instance_boundaries(self, mesh=None):
        '''Instances the relevant boundaries for boundary conditions.
        
            mesh: type of mesh created. None, "slab" or "tkd". (None or str)
        '''
        
        self.meshtype = mesh
        # Instance the relevant boundaries
        print("Instancing boundaries")

        if not self.periodic:
            self.gamma_in = GammaIn(
                self.dir_min_flow, self.dir_max_flow, self.tol
            )
            self.gamma_out = GammaOut(
                self.dir_min_flow, self.dir_max_flow, self.tol
            )
            self.gamma_air = GammaAir(
                self.dir_min_y, self.dir_max_y, self.dir_min_z, self.dir_max_z, self.tol
            )

            print("gamma_in.dir_min = ", self.gamma_in.dir_min)
            print("gamma_in.dir_max = ", self.gamma_in.dir_max)
            print("gamma_out.dir_min = ", self.gamma_out.dir_min)
            print("gamma_out.dir_max = ", self.gamma_out.dir_max)
            print("gamma_air.dir_min_y = ", self.gamma_air.dir_min_y)
            print("gamma_air.dir_max_y = ", self.gamma_air.dir_max_y)
            print("gamma_air.dir_min_z = ", self.gamma_air.dir_min_z)
            print("gamma_air.dir_max_z = ", self.gamma_air.dir_max_z)

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
            self.W_h = FunctionSpace(self.mesh, 'Lagrange', self.f_dim)
            self.V_h = VectorFunctionSpace(self.mesh, 'Lagrange', self.vecf_dim)
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
        meshtype: type of mesh. None, "sheet" or "tkd". (None or str)
        '''
        print("\n ############## (P) problem simulation ############")
        self.instance_boundaries(mesh=meshtype)
        print("boundaries instanced")
        self.instance_function_spaces()
        print("self.f_dim =", self.f_dim)
        print("self.vecf_dim =", self.vecf_dim)
        
        # Declare Dirichlet boundary conditions for (P)
        self.p_dbc = [
            DirichletBC(self.W_h, self.L*self.params['p_min']/(self.mu*self.U), self.gamma_out)
        ]

        # Assemble problem
        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        p = TrialFunction(self.W_h)
        v = TestFunction(self.W_h)
        a = inner(grad(p), grad(v))*dx
        F = Constant(0)*v*dx
        F += (self.L**2)/self.params["kappa"]*v*ds(1)

        # Solve problem

        self.p = Function(self.W_h)
        if self.solver == 'gmres':    
            solve(
                a == F, self.p, self.p_dbc
            )
        elif self.solver == "bicgstab":
            solve(
                a == F, self.p, self.p_dbc,
                solver_parameters={'linear_solver': 'bicgstab'}
            )      
        elif self.solver == 'mumps':
            solve(
                a == F, self.p, self.p_dbc,
                solver_parameters={'linear_solver': 'mumps'}
            )            
        else:
            solve(
                a == F, self.p, self.p_dbc
            )
            
        print("P problem solved")

        self.u = project(
            -1/(self.L**2)*self.params['kappa']*grad(self.p),
            self.V_h
        )
        print("u problem solved")

        self.u.rename("u", "blood velocity [um/s]")
        self.p.rename("p", "blood pressure [mmHg]")

        if save:
            print("Saving solutions")
            u_file = File(self.folder_path+'/p/u.pvd')
            u_file << self.u
            p_file = File(self.folder_path+'/p/p.pvd')
            p_file << self.p
            print("saved")

    def set_u(self, value=(0, 0, 0), meshtype='sheet', save=True):
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
            return(self.Dap_T_O2 - self.Dap_O2*c_HbX - self.Dap_CO2*c_HbY)*p_X - self.Da_O2*self.H_O2/self.Re *c_HbX
        elif X == 'CO2':
            return(self.Dap_T_CO2 - self.Dap_O2*c_HbX - self.Dap_CO2*c_HbY)*p_X - self.Da_CO2*self.H_CO2/self.Re *c_HbX
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
            return(self.S_O2_T - self.S_O2_O2*c_HbX - self.S_O2_CO2*c_HbY)*p_X - self.W_O2*c_HbX
        elif X == 'CO2':
            return(self.S_CO2_T - self.S_CO2_O2*c_HbX - self.S_CO2_CO2*c_HbY)*p_X - self.W_CO2*c_HbX
        else:
            raise ValueError('Gas species in f must be O2 or CO2.')

    def sim_t(self, hb=True, save=True, guess=None, solver=None, preconditioner=None, opt=None, stab_param=None):
        '''Solves the steady state blood-side transport (T) problem of the
        model.
        
        hb: toggle effects of hemoglobin, that is, chooses between solving
        (WT) or (WLT). (bool)
        save: saves to vtk. (bool)
        guess: starting point for Newton iterations. s^0_L in the paper.
        (FEniCS Function or None)
        '''
        
        print("\n ############## (LT) and (T) problem simulation ############")
        
        # Instance parameters
        if solver is None:
            print("No solver selected.")
        else:
            if preconditioner is None:
                print("began sim_t, with solver = ", solver, " and WITHOUT preconditioner.")
            else:
                print("began sim_t, with solver = ", solver, " and preconditioner = ", preconditioner)
        p_air_O2 = self.params['p_air_O2']
        d_ba_O2 = self.params['d_ba_O2']
        d_pla_O2 = self.params['d_pla_O2']

        p_air_CO2 = self.params['p_air_CO2']
        d_ba_CO2 = self.params['d_ba_CO2']
        d_pla_CO2 = self.params['d_pla_CO2']

        h_ba = self.params['h_ba']

        # Instance function space for the multi-field problem
        element = VectorElement('P', tetrahedron, 1, dim=4) # Era 3 previamente
        self.M_h = FunctionSpace(self.mesh, element)  

        ds = Measure('ds', domain=self.mesh, subdomain_data=self.boundaries)

        # Declare functions and test functions, and initial point for iterations

        if not guess:
            x = Function(self.M_h)
        else:
            guess.set_allow_extrapolation(True)
            x = project(guess, self.M_h)

        self.p_O2, self.p_CO2, self.c_HbO2, self.c_HbCO2 = split(x)

        v, w, eta, xi = TestFunctions(self.M_h)

        n = FacetNormal(self.mesh)
        
        # Define residuals 

        G_p_O2 = inner(grad(self.p_O2), grad(v))*dx
        G_p_O2 += -self.R_O2/self.Tau*Constant(p_air_O2)*v*ds(3)
        G_p_O2 += self.R_O2*self.p_O2*v*ds(3)
        G_p_O2 += self.Pe_O2*inner(self.p_O2*self.u, n)*v*ds(2)
        G_p_O2 += -self.Pe_O2*inner(self.p_O2*self.u, grad(v))*dx # Término convectivo
        if hb:
            G_p_O2 += self.f('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*v*dx

        G_p_CO2 = inner(grad(self.p_CO2), grad(w))*dx
        G_p_CO2 += -self.R_CO2/self.Tau*Constant(p_air_CO2)*w*ds(3)
        G_p_CO2 += self.R_CO2*self.p_CO2*w*ds(3)
        G_p_CO2 += self.Pe_O2*inner(self.p_CO2*self.u, n)*w*ds(2)
        G_p_CO2 += -self.Pe_O2*inner(self.p_CO2*self.u, grad(w))*dx # Término convectivo
        if hb:
            G_p_CO2 += self.f('CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2)*w*dx

        if hb:
            G_c_O2 = self.g('O2', self.p_O2, self.c_HbO2, self.c_HbCO2)*eta*dx
            G_c_O2 += inner(self.c_HbO2*self.u, grad(eta))*dx # Término convectivo
            G_c_O2 += -inner(self.c_HbO2*self.u, n)*eta*ds(2)

            G_c_CO2 = self.g('CO2', self.p_CO2, self.c_HbCO2, self.c_HbO2)*xi*dx
            G_c_CO2 += inner(self.c_HbCO2*self.u, grad(xi))*dx # Término convectivo
            G_c_CO2 += -inner(self.c_HbCO2*self.u, n)*xi*ds(2)

        else:
            G_c_O2 = self.c_HbO2*eta*dx
            G_c_CO2 = self.c_HbCO2*xi*dx

        G = G_p_O2 + G_p_CO2 + G_c_O2 + G_c_CO2      

        if save:
            if guess is None:
                p_O2_file = File(self.folder_path+'/t/pO2_linear.pvd')
                p_CO2_file = File(self.folder_path+'/t/pCO2_linear.pvd')
                c_HbO2_file = File(self.folder_path+'/t/cHbO2_linear.pvd')
                c_HbCO2_file = File(self.folder_path+'/t/cHbCO2_linear.pvd')
                
            else:
                p_O2_file = File(self.folder_path+'/t/pO2.pvd')
                p_CO2_file = File(self.folder_path+'/t/pCO2.pvd')
                c_HbO2_file = File(self.folder_path+'/t/cHbO2.pvd')
                c_HbCO2_file = File(self.folder_path+'/t/cHbCO2.pvd')

        # Declare Dirichlet boundary conditions for (T)
        self.t_dbc = [
            DirichletBC(
                self.M_h.sub(0), Constant(self.L*self.params['p_O2_in']/(self.mu*self.U)), self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(1), Constant(self.L*self.params['p_CO2_in']/(self.mu*self.U)),
                self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(2), Constant(1),
                self.gamma_in
            ),
            DirichletBC(
                self.M_h.sub(3), Constant(1),
                self.gamma_in
            )
        ]

        # Solve variational problem
        if solver is None:
            solve(
                G == 0, x, self.t_dbc,
                solver_parameters={"newton_solver": {
                    "relative_tolerance": 1E-8,
                    "absolute_tolerance": 1E-8#, 
#                     "linear_solver": test_solver
                }}
            )

        else:
            if preconditioner is None:
                print(f"Solving with solver = {solver} and without preconditioner.")
                solve(
                    G == 0, x, self.t_dbc,
                    solver_parameters={"newton_solver": {
                        "relative_tolerance": 1E-8,
                        "absolute_tolerance": 1E-8,
                        "linear_solver": solver
                    }}
                )

            elif preconditioner == "K":
                print("Starting Krylov solver formulation")
                
                F = Constant(0)*v*dx
                b =  inner(grad(self.p_O2),grad(v))     + inner(grad(self.p_CO2),grad(w))
                b += inner(grad(self.c_HbO2),grad(eta)) + inner(grad(self.c_HbCO2),grad(xi))
                
                # Esto no funciona porque G no es bilineal, ver alternativas
                A, bb = assemble_system(G, F, self.t_dbc)
                P, btmp = assemble_system(b, F, self.t_dbc)
                solver = KrylovSolver(krylov_method, "amg")
                solver.set_operators(A,P)
                
                solver.solve(x.vector(), bb)
            
            else:
                print(f"Solving with solver = {solver} and preconditioner = {preconditioner}.")                
                solve(
                    G == 0, x, self.t_dbc,
                    solver_parameters={"newton_solver": {
                        "relative_tolerance": 1E-8,
                        "absolute_tolerance": 1E-8,
                        "linear_solver": solver, 
                        "preconditioner": preconditioner,
                        "maximum_iterations": 100000
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