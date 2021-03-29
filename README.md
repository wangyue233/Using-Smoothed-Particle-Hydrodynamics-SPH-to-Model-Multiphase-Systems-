# ACSE-4-SPH

[Smoothed Particle Hydrodynamics](https://en.wikipedia.org/wiki/Smoothed-particle_hydrodynamics) (SPH) is a meshless
method for solving the Navier-Stokes equation, in which fluid properties are stored on Lagrangian fluid particles (i.e. on
particles which move with the fluid flow). The particles interact to generate values across the entire fluid domain through
continuous smoothing kernels. 

As the SPH method is meshless and Lagrangian, it is ideal for solving problems involving fluid flow with interfaces and free 
surfaces. This tool implements the SPH method in C++ to solve wave generation in a lock-release/dam-break problem.

### Compilation/Installation Guide

Windows, Visual studio 2017/2019:

Build a new project, add cpp files from included folder as header files, add cpp files from src files as source files.

Linux, G++ complier:

Download the whole project,  to install the code, simply run this in your favorite terminal:

make


### User instructions

SPH simulator uses two different methods(full time and half step euler method) to solve wave generation. As two methods are implemented in separate files, users are allowd to choose the method they like to realise the simulation and generate results from output folder. Besider, this simulator also enables users to decide their own simulation parameters, like: dx(space step), dt(time step), simulation time, number of output files. Output files are designed to be VTP files, which enables users to check the animation under Paraview and also do the analysis through Python. Moreover, the simulator has test files to ensure the stability for the code. Users can simply excute the command:

make runtests

to run test files. For windows system, they need to choose test file as main file to get results.

### Documentation

The full code documentation is available in html format in the ``Smooth Particles Hydrodynamics Solver Code Documentation.zip`` file. To use, unzip and drag the unzipped folder to a browser tab of your choice.

### Testing

The tool includes tests, which you can use to check its operation on your system. With the code compiled, these can be run 
with

```
python run_tests.py
make runtests
```
