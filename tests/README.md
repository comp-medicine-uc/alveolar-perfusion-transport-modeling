## Tests

- ```test01_perfusion_gas_exchange.ipynb```: Jupyter notebook tests for single-capillary slab perfusion and gas exchange results, which generate the figures shown in the article.
- ```test02_rve_DL_parallel.py```: Python script test that simulates perfusion and gas exchange in all four RVE geometries. Parallel execution via ```mpirun -n 8 python3 test02_rve_DL_parallel.py``` is recommended.
- ```./csv-results```: ```.csv``` results of Zurita & Hurtado and Brighenti et al. gas exchange simulations, used for validation plots shown in ```test01_perfusion_gas_exchange.ipynb```.
