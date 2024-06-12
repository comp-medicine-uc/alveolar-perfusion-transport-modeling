import numpy as np
import matplotlib.pyplot as plt
import pyvista
import ufl
import time
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx
from dolfinx import fem, mesh, la, plot, nls, log, io
from dolfinx import cpp as _cpp
import dolfinx.fem.petsc
import meshio
import os
import sys

sys.path.append('../')
# from src.model import PerfusionGasExchangeModel

def plot_mesh(domain):

    # Start PyVista Plotter object
    pyvista.set_jupyter_backend('client')
    pyvista.start_xvfb()
    plotter = pyvista.Plotter()

    # Create VTK mesh for visualization
        # topology = vtk topology data
        # cells = cell types
        # geometry = Nx3 array of all point coordinates

    topology, cells, geometry = plot.vtk_mesh(domain, domain.topology.dim)
    grid = pyvista.UnstructuredGrid(topology, cells, geometry)
        
    plotter.add_mesh(grid, show_edges=True)
    plotter.view_xy()
    plotter.show()

def plot_scalar_field(function_space, scalar_field):

    # Start PyVista Plotter object
    pyvista.set_jupyter_backend('client')
    pyvista.start_xvfb()
    plotter = pyvista.Plotter()
    
    # Set scalar function values
    topology, cell_types, geometry = plot.vtk_mesh(function_space)
    grid = pyvista.UnstructuredGrid(topology, cell_types, geometry)
    grid.point_data["p"] = scalar_field.x.array
    grid.set_active_scalars("p")

    plotter.add_text("field_values", position="upper_edge", font_size=14, color="black")
    plotter.add_mesh(grid, show_edges=True)
    plotter.view_xy()
    plotter.show()

def plot_vector_field(domain, vector_function_space, vector_field):

    # Start PyVista Plotter object
    pyvista.set_jupyter_backend('client')
    pyvista.start_xvfb()
    plotter = pyvista.Plotter()

    # Set vector function values
    topology1, cell_types1, geometry1 = plot.vtk_mesh(vector_function_space)
    values1 = np.zeros((geometry1.shape[0], 3), dtype=np.float64)
    values1[:, :len(vector_field)] = vector_field.x.array.real.reshape((geometry1.shape[0], len(vector_field)))
    # Create a point cloud of glyphs
    # function_grid1 = pyvista.UnstructuredGrid(topology1, cell_types1, geometry1)
    # function_grid1["u"] = values1

    # Create a pyvista-grid for the mesh
    grid = pyvista.UnstructuredGrid(*plot.vtk_mesh(domain, domain.topology.dim))
    grid.point_data["u"] = values1
    geom = pyvista.Arrow()
    glyphs = grid.glyph(orient="u", factor=0.008*10, geom=geom)

    plotter.add_mesh(grid, style="wireframe", color="k")
    plotter.add_mesh(glyphs)
    plotter.view_xy()
    plotter.show()

def mesh_to_xdmf(domain, folder, name='mesh.xdmf'):

    # Path setup
    path = os.path.join(folder, name)

    # XDMF mesh creation
    with io.XDMFFile(domain.comm, path, "w", encoding=io.XDMFFile.Encoding.ASCII) as file:
        file.write_mesh(domain)

def scalar_field_to_xdmf(domain, scalar_field, folder, name='mesh.xdmf', field_name = 'Pressure'):

    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)
        print("Folder did not exist, created.")

    with io.XDMFFile(domain.comm, os.path.join(folder, name), "w") as file:
        file.write_mesh(domain)
        # scalar_field.name = field_name
        file.write_function(scalar_field)

def vector_field_to_xdmf(domain, vector_field, folder, name='mesh.xdmf', field_name = 'Velocity'):
    
    if not os.path.exists(folder):
        os.makedirs(folder, exist_ok=True)
        print("Folder did not exist, created.")

    # XDMF mesh creation and vector field writing
    with io.XDMFFile(domain.comm, os.path.join(folder, name), "w") as file: # , encoding=io.XDMFFile.Encoding.ASCII
        file.write_mesh(domain)
        # vector_field.name = field_name
        file.write_function(vector_field)

def project(v, target_func, bcs=[]):
    # Ensure we have a mesh and attach to measure
    V = target_func.function_space
    dx = ufl.dx(V.mesh)

    # Define variational problem for projection
    w = ufl.TestFunction(V)
    Pv = ufl.TrialFunction(V)
    a = fem.form(ufl.inner(Pv, w) * dx)
    L = fem.form(ufl.inner(v, w) * dx)

    # Assemble linear system
    A = fem.petsc.assemble_matrix(a, bcs)
    A.assemble()
    b = fem.petsc.assemble_vector(L)
    fem.petsc.apply_lifting(b, [a], [bcs])
    b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE)
    fem.petsc.set_bc(b, bcs)

    # Solve linear system
    solver = PETSc.KSP().create(A.getComm())
    solver.setOperators(A)
    solver.solve(b, target_func.vector)

def L2_norm(comm, v):
    """Compute the L2(Ω)-norm of v"""
    return np.sqrt(comm.allreduce(fem.assemble_scalar(fem.form(ufl.inner(v, v) * ufl.dx)), op=MPI.SUM))

def domain_average(domain, v):
    # Compute the average of a function over the domain
    vol = domain.comm.allreduce(fem.assemble_scalar(fem.form(fem.Constant(domain, dolfinx.default_real_type(1.0)) * ufl.dx)), op=MPI.SUM)
    return (1 / vol) * domain.comm.allreduce(fem.assemble_scalar(fem.form(v * ufl.dx)), op=MPI.SUM)

def calculate_area(domain, ds, tag):
    return domain.comm.allreduce(fem.assemble_scalar(fem.form(fem.Constant(domain, dolfinx.default_real_type(1.0)) * ds(tag))), op=MPI.SUM)

def calculate_flux(domain, ds, tag, f, v, n):
    return domain.comm.allreduce(fem.assemble_scalar(fem.form(ufl.inner(v, n) * f * ds(tag))), op=MPI.SUM)

def calculate_integral_over_domain(domain, dx, f):
    return domain.comm.allreduce(fem.assemble_scalar(fem.form(f * dx)), op=MPI.SUM)

def calculate_integral_over_surface(domain, ds, tags, f):
    result = 0
    for tag in tags:
        result += domain.comm.allreduce(fem.assemble_scalar(fem.form(f * ds(tag))), op=MPI.SUM)
    return result

def diagnose_mesh(domain, dx):
    _, _, geometry = plot.vtk_mesh(domain, domain.topology.dim)
        
    max_dims = np.around([max(geometry[:,0]), max(geometry[:,1]), max(geometry[:,2])], 5)
    min_dims = np.around([min(geometry[:,0]), min(geometry[:,1]), min(geometry[:,2])], 5)

    print(f"max_dims = {np.round(max_dims, 5)}")
    print(f"min_dims = {np.round(min_dims, 5)}")
    
    tissue_volume = calculate_integral_over_domain(domain, dx, 1)
    all_volume = np.prod(max_dims - min_dims)
    air_volume = all_volume - tissue_volume

    porosity = air_volume/all_volume

    print(f"all volume = {np.round(all_volume, 5)}")
    print(f"tissue volume = {np.round(tissue_volume, 5)}")
    print(f"air volume = {np.round(air_volume, 5)}")
    print(f"porosity = {np.round(porosity, 5)}")

    return {"max_dims": max_dims, "min_dims": min_dims, "tissue_volume": tissue_volume, 
            "all_volume": all_volume, "air_volume": air_volume, "porosity": porosity}

def read_diagnosis(meshname, rank):

    metrics = []

    with open("../tests/csv-results/geometry_data.csv", "r") as file:
        for row in file:
            items = row.split(",")
            if items[0] in meshname:
                metrics.append([float(i) for i in items[1].strip("]").strip("[").strip().replace("  ", " ").replace(" ", ",").split(",")])
                metrics.append([float(i) for i in items[2].strip("]").strip("[").strip().replace("  ", " ").replace(" ", ",").split(",")])
                metrics.append(items[3])
                metrics.append(items[4])
                metrics.append(items[5])
                metrics.append(items[6])
        
    if rank == 0:
        print(f"###### Geometry metrics ######")
        print(f"max_dims = {metrics[0]}")
        print(f"min_dims = {metrics[1]}")
        print(f"all_volume = {metrics[2]}")
        print(f"tissue_volume = {metrics[3]}")
        print(f"air_volume = {metrics[4]}")
        print(f"porosity = {metrics[5]}")
    return metrics
               

def plot_over_lines(domain, field_dict: dict, y=0, z=0):
    """
    Make plots for curves at slab points. Adaptation of
    https://jsdokken.com/dolfinx-tutorial/chapter1/membrane_code.html#making-curve-plots-throughout-the-domain
    """
    _, _, geometry = plot.vtk_mesh(domain, domain.topology.dim)
        
    max_dims = np.around([max(geometry[:,0]), max(geometry[:,1]), max(geometry[:,2])], 5)
    min_dims = np.around([min(geometry[:,0]), min(geometry[:,1]), min(geometry[:,2])], 5)

    tol = 0.001
    npoints = 1000
    
    x = np.linspace(min_dims[0] + tol, max_dims[0] - tol, npoints)
    _y = y * np.ones_like(x)
    _z = z * np.ones_like(x)

    points = np.zeros((3, npoints))
    points[0] = x
    points[1] = _y
    points[2] = _z

    bb_tree = dolfinx.geometry.bb_tree(domain, domain.topology.dim)

    cells = []
    points_on_proc = []

    # Find cells whose bounding-box collide with the the points
    cell_candidates = dolfinx.geometry.compute_collisions_points(bb_tree, points.T)
    # Choose one of the cells that contains the point
    colliding_cells = dolfinx.geometry.compute_colliding_cells(domain, cell_candidates, points.T)
    for i, point in enumerate(points.T):
        if len(colliding_cells.links(i)) > 0:
            points_on_proc.append(point)
            cells.append(colliding_cells.links(i)[0])

    points_on_proc = np.array(points_on_proc, dtype=np.float64)

    all_func_names = []
    all_func_values = []
    # print(f"Evaluating...")
    for func_name, func in field_dict.items():
        # print(f"... current function: {func_name}")
        func_values = func.eval(points_on_proc, cells)
        # print("evaluated.")
        all_func_names.append(func_name)
        all_func_values.append(func_values)

    return all_func_names, all_func_values, x