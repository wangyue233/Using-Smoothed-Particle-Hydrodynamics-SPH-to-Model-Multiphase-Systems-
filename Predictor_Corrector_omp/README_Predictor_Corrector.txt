Final and improved version of the solver

This version includes the following additions/optimisations to the SPH solver:

a) Omp parallel implementation
b) Predictor-Corrector method to increase stability
c) Use of neighbor stencil to avoid applying symmetric interactions
d) Skipping the calculations of boundary particles (except pressure update)
e) Proper application of mirror boundary conditions
f) Precomputed application of Weights function to favor computation speed
g) General improvements in linting and code structure

Computation timings for varying dx and 5000 iterations and constant h factor = 1.3:

dx	dt	time (sec)
0.2	1.3e-3	12
0.1	6.5e-4	35
0.05	3.2e-4	150
0.02	1.3e-4	280

CPU: i5 8 logical cores processor
