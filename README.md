# Computational simulations of an alveolar perfusion and gas transport mathematical model

[Bastián Herrera](https://github.com/bnherrerac), _[@comp-medicine-uc](https://github.com/comp-medicine-uc)_

FEniCS-based computational simulation of alveolar perfusion and gas transport on 3D alveolar meshes.

## Abstract

(Add paper abstract here)

## Directories

- `manuscript`: Files pertaining to the article related to this project.
- `presentations`: Files pertaining to presentations of this project.
- `raw-and-results-data`: Folders pertaining to each different tests, which contain raw images + outputs.
- `src`: Source files.
- `tests`: Main files that implement examples and tests.

## Dependencies

Project tested in MATLAB R2021a Update 6, Python 3.

- [`iso2mesh`](http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Home) 2018.1.9.6 Lion's Head

Used libraries:

- [`FEniCS`](https://fenicsproject.org/) 2019.1.0
- [`TetGen bindigs for Python`](https://tetgen.pyvista.org) 0.6.1
- [`Trimesh`](https://trimsh.org/trimesh.html) 3.15.1
- [`PyVista`](https://docs.pyvista.org) 0.34.0
- [`PyMeshFix`](https://pymeshfix.pyvista.org) 0.16.1, a Python/Cython wrapper for [`MeshFix`](https://github.com/MarcoAttene/MeshFix-V2.1).
- [`numpy`](https://numpy.org) 1.23.4
- [`scipy`](https://scipy.org) 1.9.2
- [`matplotlib`](https://matplotlib.org) 3.6.1
- `os`, `sys`