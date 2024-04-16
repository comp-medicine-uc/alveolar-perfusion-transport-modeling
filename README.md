# A nonlinear convection-diffusion model of microscale lung capillary perfusion and gas exchange

[Bastián Herrera](https://github.com/bnherrerac) & [Daniel Hurtado](https://github.com/dehurtado), _[@comp-medicine-uc](https://github.com/comp-medicine-uc)_

Computational simulations of microscale pulmonary perfusion and O$`_2`$-CO$`_2`$ gas transport on 3D alveolar meshes reconstructed from $`\mu`$-CT murine lung images.

## Abstract

Pulmonary capillary perfusion and gas exchange are the fundamental physical processes that enable respiration to occur at the microscale. Present-day computational simulations of these phenomena are often based on low-dimensional mathematical models on idealized alveolar geometries, where the chemical reactions between O$`_2`$-CO$`_2`$ and haemoglobin are simplified. However, these models fail to capture the complex chemical reactions that take place in pulmonary capillary blood and overlook the effect of the porous alveolar structure. In this work, we develop a coupled perfusion and gas exchange model that reflects gas and haemoglobin dynamics in pulmonary capillaries. To this end, we formulate a system of two coupled nonlinear convection-diffusion boundary value problems. We derive the equations from continuum balance laws, incorporating an experimentally validated relationship between gas partial pressures and haemoglobin saturations through a fully coupled Hill equation-like model. We numerically solve this boundary value problem in Python using the finite element method in a simple slab capillary domain, considering a physiological range of blood flow velocity fields and blood pH values. Moreover, we perform numerical simulations in a 3D alveolar domain reconstructed from $`\mu`$-CT rat lung images, which was morphologically altered to simulate various levels of emphysema. Numerical perfusion experiments agree with expected blood pressure drops and velocity fields in the lung capillary domain. Additionally, numerical gas exchange simulations reconstruct physiological dynamics of gas partial pressure and haemoglobin saturation in the capillary domain. Furthermore, these simulations replicate the Bohr effect and provide a framework to determine the whole-lung diffusing capacity for O$`_2`$ and CO$`_2`$, which diminishes as emphysema severity increases. We envision that this model broadens potential applicability of computational lung simulations in clinical settings, particularly in exercise and pathological conditions that affect perfusion dynamics and the overall gas exchange function of the lungs.

## Directories

- `manuscript`: Latest version of the article.
- `raw-data`: Raw images and sources for pulmonary RVE simulations.}
- `results-data`: Numerical simulation results.
- `src`: Python and MATLAB source files.
- `tests`: Main files that implement examples and tests.

## Dependencies

Source code is written and runs in Python 3.12 and in MATLAB R2023a. The following libraries are employed:

- [`iso2mesh`](http://iso2mesh.sourceforge.net/cgi-bin/index.cgi?Home) for MATLAB, 2018.1.9.6 Lion's Head
- [`numpy`](https://numpy.org) 1.26.0
- [`scipy`](https://scipy.org) 1.13.0
- [`matplotlib`](https://matplotlib.org) 3.8.4
- [`FEniCSx`](https://fenicsproject.org/) 0.7.2
- [`TetGen bindigs for Python`](https://tetgen.pyvista.org) 0.6.1
- [`Trimesh`](https://trimsh.org/trimesh.html) 4.3.0
- [`PyVista`](https://docs.pyvista.org) 0.43.5
- [`mpi4py`](https://mpi4py.readthedocs.io/en/stable/index.html) 3.1.3
- [`PyMeshFix`](https://pymeshfix.pyvista.org) 0.16.3, a Python/Cython wrapper for [`MeshFix`](https://github.com/MarcoAttene/MeshFix-V2.1).
- [`meshio`](https://pypi.org/project/meshio/) 5.3.5
- [`petsc4py`](https://petsc.org/release/petsc4py/) 3.15.1
- `os`, `sys`
