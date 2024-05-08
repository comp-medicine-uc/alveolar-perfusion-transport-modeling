import numpy as np
import matplotlib.pyplot as plt
import pyvista
import ufl
import time
from mpi4py import MPI
from petsc4py import PETSc
import dolfinx
from dolfinx import fem, mesh, plot, nls, log, io
from dolfinx import cpp as _cpp
import meshio
import os
import sys
import csv

sys.path.append('../')
sys.path.append('./')

from src.utils import *
from src.model import *

with open("./csv-results/geometry_data.csv", "w") as file:
    writer = csv.writer(file)
    for root, dirs, files in os.walk("../"):
        for file in files:
            if file.endswith(".xdmf") and "results-data" not in root:
                path = os.path.join(root, file)
                print(path)

                with io.XDMFFile(MPI.COMM_WORLD, path, 'r') as xdmf:
                    domain = xdmf.read_mesh(name='Grid')

                    metadata = {"quadrature_degree": 4}
                    diagnosis = diagnose_mesh(domain, ufl.Measure("dx", domain=domain, metadata=metadata))
                    print("------------------------------------------------------")

                    row = [file] + list(diagnosis.values())

                    writer.writerow(row)
