## Raw data

This folder contains raw micro-CT rat lung images, and RVE meshes generated using ```mesh_generation.m``` and ```mesh_processing.py```.

- ```./binary_slices```: binary micro-CT rat lung slices.
- ```./rve_225_emphysema_nodes_elems```: node and element arrays for all 4 RVE meshes (normal and emphysema 1-2-3).
- ```./rve_225_emphysema_control```: control, non-eroded RVE mesh reconstructed from the binary slices.
- ```./rve_225_emphysema_1```: emphysematous RVE mesh, eroded with a small disk.
- ```./rve_225_emphysema_2```: emphysematous RVE mesh, eroded with a medium disk.
- ```./rve_225_emphysema_3```: emphysematous RVE mesh, eroded with a large disk.

The full meshing and processing algorithms are detailed in ```mesh_generation.m``` and ```mesh_processing.py```.
