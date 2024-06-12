import numpy as np
import matplotlib.pyplot as plt
import pyvista
import ufl
import time
from mpi4py import MPI
from petsc4py import PETSc
from dolfinx import fem, mesh, plot, nls, log, io, la
import dolfinx.fem.petsc
import dolfinx.nls.petsc
import meshio
import os
import sys

sys.path.append('../')

from src.utils import *
from src.utils import project
from src.db_sats import *
from src.parameters import instance_params

class PerfusionGasExchangeModel():

    def __init__(self, mesh_path, results_path, exp_path, params, ksp_type='cg', pc_type='lu', pc_factor_mat_solver_type='mumps'):

        # -pc_factor_mat_solver_type - petsc, superlu, superlu_dist, mumps, cusparse

        self.mesh_path = mesh_path
        self.results_path = os.path.join(results_path, exp_path)
        self.params = params
        self.ksp_type = ksp_type
        self.pc_type = pc_type
        self.pc_factor_mat_solver_type = pc_factor_mat_solver_type

    def Setup(self, domain, atol=1E-10, max_dims=[0,0,0], min_dims=[0,0,0], imported=True, infinite=False):

        metadata = {"quadrature_degree": 4}

        comm = MPI.COMM_WORLD
        self_comm = MPI.COMM_SELF
        comm_size = comm.Get_size()
        rank = comm.Get_rank()

        if imported:
            diagnosis = read_diagnosis(self.mesh_path, MPI.COMM_WORLD.Get_rank())
            self.max_dims = diagnosis[0]
            self.min_dims = diagnosis[1]
            self.porosity = float(diagnosis[5])
        else:
            self.max_dims = max_dims
            self.min_dims = min_dims
            self.porosity = 0

        assert np.amax(np.concatenate((self.max_dims, self.min_dims))) != 0

        # if MPI.COMM_WORLD.Get_rank() == 0: 
        #     print(f"max_dims = {self.max_dims}")
        #     print(f"min_dims = {self.min_dims}")

        self.atol = atol
        self.imported = imported

        def inlet(x):
            return np.isclose(x[0], self.min_dims[0], atol=self.atol)
        def outlet(x):
            return np.isclose(x[0], self.max_dims[0], atol=self.atol)
        def all(x):
            return np.full(x.shape[1], True, dtype=bool)
        def both(x):
            return np.logical_or(inlet(x), outlet(x))
        # def sides(x):
        #     return np.logical_or(np.isclose(x[2], self.min_dims[2], atol=self.atol), np.isclose(x[2], self.max_dims[2], atol=self.atol))

        self.fdim = domain.topology.dim - 1

        inlet_facets = mesh.locate_entities_boundary(domain, self.fdim, inlet) # IN
        outlet_facets = mesh.locate_entities_boundary(domain, self.fdim, outlet)   # OUT
        all_facets = mesh.locate_entities_boundary(domain, self.fdim, all)
        both_facets = mesh.locate_entities_boundary(domain, self.fdim, both)
        air_facets = np.setdiff1d(all_facets, both_facets)

        marked_facets = np.hstack([inlet_facets, outlet_facets, air_facets])
        marked_values = np.hstack([np.full_like(inlet_facets, 1), 
                                np.full_like(outlet_facets, 2),
                                np.full_like(air_facets, 3)])
        
        # Face numbers

        tfs = comm.gather(all_facets, root=0)
        ifs = comm.gather(inlet_facets, root=0)
        ofs = comm.gather(outlet_facets, root=0)
        afs = comm.gather(air_facets, root=0)

        if rank == 0:

            tfn = np.concatenate(tuple(tfs[i] for i in range(comm_size)), axis=None).shape[0]
            ifn = np.concatenate(tuple(ifs[i] for i in range(comm_size)), axis=None).shape[0]
            ofn = np.concatenate(tuple(ofs[i] for i in range(comm_size)), axis=None).shape[0]
            afn = np.concatenate(tuple(afs[i] for i in range(comm_size)), axis=None).shape[0]

            print(f"Porosity = {np.round(self.porosity, 3)}")
            print(f"Total boundary face number = {tfn}")
            print(f"Inlet face number = {ifn}")
            print(f"Outlet face number = {ofn}")
            print(f"Air face number = {afn}")

            # Assert that inlet and outlet are well-defined.
            assert ifn != 0
            assert ofn != 0

        sorted_facets = np.argsort(marked_facets)
        self.facet_tag = mesh.meshtags(domain, self.fdim, marked_facets[sorted_facets], marked_values[sorted_facets])

        self.ds = ufl.Measure('ds', domain=domain, subdomain_data=self.facet_tag, metadata=metadata)
        self.dx = ufl.Measure("dx", domain=domain, metadata=metadata)
        self.n = ufl.FacetNormal(domain)

    def parameter_setup(self):

        perfusion_params, dash_params, transport_params = instance_params()
        self.p_params = perfusion_params
        self.dash_params = dash_params
        self.t_params = transport_params

    def Perfusion(self, domain, plot=True, save=True, verbose=True):
        
        ### Problem setup

        # Taylor-Hood spaces
        V = fem.FunctionSpace(domain, ("CG", 2))
        W = fem.VectorFunctionSpace(domain, ("CG", 1))

        p_params = self.p_params

        metadata = {"quadrature_degree": 4}

        ds = ufl.Measure('ds', domain=domain, subdomain_data=self.facet_tag, metadata=metadata)
        dx = ufl.Measure("dx", domain=domain, metadata=metadata)
        
        # if self.imported:
        outlet_dofs = fem.locate_dofs_topological(V, self.facet_tag.dim, self.facet_tag.find(2))
        pressure_bc = fem.dirichletbc(PETSc.ScalarType(p_params["pmin"]), outlet_dofs, V)

        # Trial and test functions
        p = ufl.TrialFunction(V)
        v = ufl.TestFunction(V)

        # Linear and bilinear form assembly
        a = ufl.dot(ufl.grad(p), ufl.grad(v)) * dx
        L = PETSc.ScalarType((1/p_params["kappa"])*p_params["mu"]*p_params["uin"]) * v * ds(1)

        # Problem instancing
        problem = dolfinx.fem.petsc.LinearProblem(a, L, bcs=[pressure_bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
        if MPI.COMM_WORLD.Get_rank() == 0: print("-------   Perfusion problem instanced.  -------")

        ### Solving

        # Pressure
        if verbose:
            log.set_log_level(log.LogLevel.INFO)
        else:
            log.set_log_level(log.LogLevel.ERROR)
        ph = problem.solve()
        log.set_log_level(log.LogLevel.ERROR)
        if MPI.COMM_WORLD.Get_rank() == 0: print("Pressure field found.")
        if plot: plot_scalar_field(V, ph)

        # Velocity
        u = fem.Function(W)
        project(-1/p_params['mu']*p_params['kappa']*ufl.grad(ph),
                u)
        
        if plot: plot_vector_field(domain, W, u)

        ## Saving

        if save:
            # Interpolation needed because mesh has degree = 1.
            p_save = fem.Function(fem.FunctionSpace(domain, ("CG", 1)))
            p_save.interpolate(ph)
            
            scalar_field_to_xdmf(domain, p_save, 
                                folder = os.path.join(self.results_path, "p_new"), 
                                name = "p.xdmf")

            vector_field_to_xdmf(domain, u, 
                                folder = os.path.join(self.results_path, "u_new"), 
                                name = "u.xdmf")
        
        # Assign fields to model

        self.p = ph
        self.u = u

        return ph, u
    
    def GasExchange(self, domain, guess=None, save=True, plot=True, p_val = 1, postprocess=True, plot_lines=True, verbose=True):

        # Two Lagrange P1 spaces are used.        

        element_type = "CG"
        self.element_type = element_type
        # print("Finite element is " + element_type)

        P1 = ufl.FiniteElement(element_type, domain.ufl_cell(), 1)
        ME = ufl.MixedElement([P1, P1])
        V = fem.FunctionSpace(domain, ME)
        
        # Unknown fields
        p_XY = fem.Function(V)
        p_O2, p_CO2 = ufl.split(p_XY)

        if guess is not None:
            # p_XY = project(guess, p_XY)
            p_XY.split()[0].interpolate(guess.sub(0))
            p_XY.split()[1].interpolate(guess.sub(1))
            # print(domain_average(domain, p_XY.split()[0]))
            # print(domain_average(domain, p_XY.split()[1]))
        
        # Test functions
        v, w = ufl.TestFunctions(V)
        
        metadata = {'quadrature_degree': 4}
        ds = ufl.Measure('ds', domain=domain, subdomain_data=self.facet_tag, metadata=metadata)
        dx = ufl.Measure("dx", domain=domain, metadata=metadata)

        # Parámetros definidos en Dash & Bassingthwaighte (2010, 2016).
        W_bl = self.dash_params["W_bl"]             # adimensional, con Hct = 0.45
        beta_O2 = self.dash_params["beta_O2"]       # mmol/(um3*mmHg)
        beta_CO2 = self.dash_params["beta_CO2"]     # mmol/(um3*mmHg)
        Hb_bl = self.dash_params["Hb_bl"]           # mmol/um3
        Hct = self.dash_params["Hct"]               # adimensional (fracción volumétrica de rbcs en sangre)
        W_rbc = self.dash_params["W_rbc"]           # adimensional
        W_pl = self.dash_params["W_pl"]             # adimensional
        K1 = self.dash_params["K1"]                 # mmol/um3 (constante de hidratación CO2 + H2O = HCO3 + H)
        pH_rbc = self.dash_params["pH_rbc"]         # log(mmol/um3)
        pH_pla = self.dash_params["pH_pla"]
        c_H = 10**(-pH_rbc)*1E-12                   # mmol/um3
        R_rbc = 10**(1.357-0.205*pH_pla)
        bic = ((1-Hct)*W_pl+Hct*W_rbc*R_rbc)*(K1*beta_CO2/c_H)
    
        # Saturaciones
        S_HbO2, S_HbCO2 = S_HbXY(p_O2, p_CO2, self.dash_params)
        # print(f"domain average S_HbO2 = {domain_average(domain, S_HbO2)}")
        # print(f"domain average S_HbCO2 = {domain_average(domain, S_HbCO2)}")        

        # beta_O2 en mmol/(L*mmHg), p_O2 en mmHg, c_O2 en mmol/L = mM

        # Concentraciones
        c_O2 = ufl.variable(W_bl*beta_O2*p_O2)                      # W_bl*beta_O2*p_O2 + 4*Hb_bl*S_HbO2
        c_CO2 = ufl.variable(W_bl*beta_CO2*p_CO2 + bic*p_CO2)       # W_bl*beta_CO2*p_CO2 + bic*p_CO2 + 4*Hb_bl*S_HbCO2

        if MPI.COMM_WORLD.Get_rank() == 0: print(f"------- Gas exchange problem instanced. ------- lambda = {np.round(p_val, 3)}")

        if guess is not None:
            # print(fr"Starting nonlinear problem with previous guess")
            c_O2 += p_val*4*Hb_bl*S_HbO2
            c_CO2 += p_val*4*Hb_bl*S_HbCO2

        # Defining F
        n = ufl.FacetNormal(domain)

        d_pla_O2 = self.t_params["d_pla_O2"]                # um2/s
        d_pla_CO2 = self.t_params["d_pla_CO2"]              # um2/s
        d_ba_O2 = self.t_params["d_ba_O2"]                  # um2/s
        d_ba_CO2 = self.t_params["d_ba_CO2"]                # um2/s
        h_ba = self.t_params["h_ba"]                        # um
        p_O2_air = self.t_params["p_O2_air"]                # mmHg
        p_CO2_air = self.t_params["p_CO2_air"]              # mmHg
        p_O2_in = self.t_params["p_O2_in"]                  # mmHg
        p_CO2_in = self.t_params["p_CO2_in"]                # mmHg        
        
        # Oxygen
        F_O2 =  d_pla_O2 * ufl.inner(ufl.grad(c_O2), ufl.grad(v)) * dx              # Dif 1
        F_O2 += (d_ba_O2 / h_ba) * v * (beta_O2*p_O2) * ds(3)                       # Dif 2 (variable)
        F_O2 -= (d_ba_O2 / h_ba) * v * beta_O2 * p_O2_air * ds(3)                   # Dif 3 (constante)
        F_O2 += ufl.inner(self.u * c_O2, n) * v * ds(2)                             # Adv 1
        F_O2 -= ufl.inner(self.u * c_O2, ufl.grad(v)) * dx                          # Adv 2

        # Carbon dioxide
        F_CO2 =  d_pla_CO2 * ufl.inner(ufl.grad(c_CO2), ufl.grad(w)) * dx           # Dif 1
        F_CO2 += (d_ba_CO2 / h_ba) * w * (beta_CO2*p_CO2) * ds(3)                   # Dif 2 (variable)
        F_CO2 -= (d_ba_CO2 / h_ba) * w * beta_CO2 * p_CO2_air * ds(3)               # Dif 3 (constante)
        F_CO2 += ufl.inner(self.u * c_CO2, n) * w * ds(2)                           # Adv 1
        F_CO2 -= ufl.inner(self.u * c_CO2, ufl.grad(w)) * dx                        # Adv 2

        # Complete functional
        F = F_O2 + F_CO2

        # Dirichlet boundary condition (pressure prescribed at inlet)
        inlet_dofs_O2 = fem.locate_dofs_topological(V.sub(0), self.facet_tag.dim, self.facet_tag.find(1))
        inlet_dofs_CO2 = fem.locate_dofs_topological(V.sub(1), self.facet_tag.dim, self.facet_tag.find(1))
        bc_O2 = fem.dirichletbc(dolfinx.default_scalar_type(p_O2_in), inlet_dofs_O2, V.sub(0))
        bc_CO2 = fem.dirichletbc(dolfinx.default_scalar_type(p_CO2_in), inlet_dofs_CO2, V.sub(1))
        bcs = [bc_O2, bc_CO2]


        # # Nonlinear problem - KSP
        problem = dolfinx.fem.petsc.NonlinearProblem(F, p_XY, bcs)
        solver = dolfinx.nls.petsc.NewtonSolver(domain.comm, problem)
        solver.atol = 1E-8
        solver.rtol = 1E-8
        solver.convergence_criterion = "incremental"
        ksp = solver.krylov_solver
        opts = PETSc.Options()
        option_prefix = ksp.getOptionsPrefix()
        try:
            opts[f"{option_prefix}ksp_type"] = self.ksp_type
            opts[f"{option_prefix}pc_type"] = self.pc_type
            opts[f"{option_prefix}pc_factor_mat_solver_type"] = self.pc_factor_mat_solver_type
        except:
            raise ValueError('Option not allowed.')
        ksp.setFromOptions()

        # Solve problem
        if verbose:
            log.set_log_level(log.LogLevel.INFO)
        else:
            log.set_log_level(log.LogLevel.ERROR)
        num_its, converged = solver.solve(p_XY)
        assert(converged)
        if MPI.COMM_WORLD.Get_rank() == 0: print(f"O2 and CO2 partial pressures found in {num_its} iterations.")


        ## Postprocessing
        log.set_log_level(log.LogLevel.ERROR)
        _p_O2, _p_CO2 = p_XY.split()
        
        # if plot:
        #     plot_scalar_field(V.sub(0).collapse()[0], _p_O2)
        #     plot_scalar_field(V.sub(1).collapse()[0], _p_CO2)

        V_mesh = fem.functionspace(domain, ('Lagrange', 1))

        p_O2_solution = fem.Function(V_mesh, name="p_O2")
        p_O2_solution.interpolate(fem.Expression(_p_O2, V_mesh.element.interpolation_points()))

        p_CO2_solution = fem.Function(V_mesh, name="p_CO2")
        p_CO2_solution.interpolate(fem.Expression(_p_CO2, V_mesh.element.interpolation_points()))

        if save:
            scalar_field_to_xdmf(domain, p_O2_solution, 
                        folder = os.path.join(self.results_path, "gas_exchange/"), 
                        name = "p_O2.xdmf")
            
            scalar_field_to_xdmf(domain, p_CO2_solution, 
                        folder = os.path.join(self.results_path, "gas_exchange/"), 
                        name = "p_CO2.xdmf")
            
        if postprocess:
            if MPI.COMM_WORLD.Get_rank() == 0: print(f"Starting postprocessing.")
            t1 = time.time()

            # S_HbO2
            S_HbO2_solution = fem.Function(V_mesh, name="S_HbO2")
            S_HbO2_solution.interpolate(fem.Expression(S_HbO2, V_mesh.element.interpolation_points()))
            scalar_field_to_xdmf(domain, S_HbO2_solution, 
                        folder = os.path.join(self.results_path, "gas_exchange/"), 
                        name = "S_HbO2.xdmf")

            # S_HbCO2
            S_HbCO2_solution = fem.Function(V_mesh, name="S_HbCO2")
            S_HbCO2_solution.interpolate(fem.Expression(S_HbCO2, V_mesh.element.interpolation_points()))
            scalar_field_to_xdmf(domain, S_HbCO2_solution, 
                        folder = os.path.join(self.results_path, "gas_exchange/"), 
                        name = "S_HbCO2.xdmf")
            
            # # grad(p_O2)
            # grad_p_O2_out = fem.Function(fem.VectorFunctionSpace(domain, ('CG',1)), name="grad(p_O2)") 
            # project(ufl.grad(p_O2), grad_p_O2_out)

            # grad_p_O2_xdmf = io.XDMFFile(domain.comm, os.path.join(self.results_path, "gas_exchange/grad_p_O2.xdmf"), "w")
            # grad_p_O2_xdmf.write_mesh(domain)
            # grad_p_O2_xdmf.write_function(grad_p_O2_out, 0.)

            # # grad(p_CO2)
            # grad_p_CO2_out = fem.Function(fem.VectorFunctionSpace(domain, ('CG',1)), name="grad(p_CO2)") 
            # project(ufl.grad(p_CO2), grad_p_CO2_out)

            # grad_p_CO2_xdmf = io.XDMFFile(domain.comm, os.path.join(self.results_path, "gas_exchange/grad_p_CO2.xdmf"), "w")
            # grad_p_CO2_xdmf.write_mesh(domain)
            # grad_p_CO2_xdmf.write_function(grad_p_CO2_out, 0.)

            # # j_O2 = -d_pla_O2 * grad(c_O2) (diffusive flux)
            # j_O2_out = fem.Function(fem.VectorFunctionSpace(domain, ('CG',1)), name="j_O2") 
            # project(-d_pla_O2 * ufl.grad(c_O2), j_O2_out)

            # j_O2_xdmf = io.XDMFFile(domain.comm, os.path.join(self.results_path, "gas_exchange/j_O2.xdmf"), "w")
            # j_O2_xdmf.write_mesh(domain)
            # j_O2_xdmf.write_function(j_O2_out, 0.)

            # # j_CO2 = -d_pla_CO2 * grad(c_CO2) (diffusive flux)
            # j_CO2_out = fem.Function(fem.VectorFunctionSpace(domain, ('CG',1)), name="j_CO2") 
            # project(-d_pla_CO2 * ufl.grad(c_CO2), j_CO2_out)

            # j_CO2_xdmf = io.XDMFFile(domain.comm, os.path.join(self.results_path, "gas_exchange/j_CO2.xdmf"), "w")
            # j_CO2_xdmf.write_mesh(domain)
            # j_CO2_xdmf.write_function(j_CO2_out, 0.) 

            # c_O2
            c_O2_solution = fem.Function(V_mesh, name="c_O2")
            c_O2_solution.interpolate(fem.Expression(c_O2, V_mesh.element.interpolation_points()))
            scalar_field_to_xdmf(domain, c_O2_solution, 
                        folder = os.path.join(self.results_path, "gas_exchange/"), 
                        name = "c_O2.xdmf")

            # c_CO2
            c_CO2_solution = fem.Function(V_mesh, name="c_CO2")
            c_CO2_solution.interpolate(fem.Expression(c_CO2, V_mesh.element.interpolation_points()))
            scalar_field_to_xdmf(domain, c_CO2_solution, 
                        folder = os.path.join(self.results_path, "gas_exchange/"), 
                        name = "c_CO2.xdmf")

            # bicarbonate CO2 content
            bic_solution = fem.Function(V_mesh, name="bic")
            bic_solution.interpolate(fem.Expression(bic*p_CO2, V_mesh.element.interpolation_points()))
            scalar_field_to_xdmf(domain, bic_solution, 
                        folder = os.path.join(self.results_path, "gas_exchange/"), 
                        name = "bic.xdmf")


            if MPI.COMM_WORLD.Get_rank() == 0: print(f"Finished postprocessing at t = {np.round(time.time() - t1, 4)} s.")   

            self.S_HbO2 = S_HbO2_solution
            self.S_HbCO2 = S_HbCO2_solution
            self.c_O2 = c_O2_solution
            self.c_CO2 = c_CO2_solution
            self.bic = bic_solution
            self.p_O2 = p_O2_solution
            self.p_CO2 = p_CO2_solution

        if plot_lines:
                
            _, _, geometry = dolfinx.plot.vtk_mesh(domain, domain.topology.dim)
            max_dims = np.around([max(geometry[:,0]), max(geometry[:,1]), max(geometry[:,2])], 5)
            min_dims = np.around([min(geometry[:,0]), min(geometry[:,1]), min(geometry[:,2])], 5)
            zs = np.linspace((max_dims[2]+min_dims[2])/2, max_dims[2], 3)
            y_fixed = (max_dims[1]+min_dims[1])/2

            func_dict = {
                "S_HbO2": S_HbO2_solution,
                         "S_HbCO2": S_HbCO2_solution,
                         "c_O2": c_O2_solution, 
                         "c_CO2": c_CO2_solution,
                         "bic": bic_solution,
                         "p_O2": p_O2_solution,
                         "p_CO2": p_CO2_solution}

            results_dict = {key: {} for key in func_dict.keys()} # {bic: {4: results_4, 6: results_6, 8: results_8}}

            for z in zs:
                all_func_names_x, all_func_values_x, x = plot_over_lines(domain, func_dict, y=y_fixed, z=z)
                for i in range(len(all_func_values_x)):
                    results_dict[all_func_names_x[i]][f"{z}"] = all_func_values_x[i]

            self.results_dict = results_dict
            self.func_dict = func_dict
            self.x = x
            self.zs = zs
            
        return p_XY