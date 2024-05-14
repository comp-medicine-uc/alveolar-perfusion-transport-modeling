## Source files

This folder contains the model source files and utility functions.

- ```db_sats.py```: Python functions for Dash & Bassingthwaighte expressions of $S_{HbO_2}$ and $S_{HbCO}_2$ as a function of both gas partial pressures.
- ```mesh_generation.m```: MATLAB image processing algorithm for surface mesh generation from binary slices.
- ```img_merge.m```: MATLAB script that joins slices into a 3D array.
- ```mesh_processing.py```: Python script that cleans and repairs surface meshes, tetrahedralizes them into volume meshes, and outputs a ```h5``` mesh.
- ```mesh_diagnosis.py```: Python script that returns geometry info from a mesh.
- ```model.py```: Python model implementation.
- ```parameters.py```: Model parameter instancing.
- ```utils.py```: Python utils for finite element analysis and calculations performed by the model.
- ```diffusing_capacity_utils.py```: Python utils for diffusing capacity calculations from microscale gas exchange results.
- ```./fonts```: TrueType font files used for figure generation.