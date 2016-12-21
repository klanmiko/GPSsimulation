GPS Simulator
==============
GPS.m is a MATLAB script that defines a function GPS(target)
target:
   vector of [x y z] coordinates to test the accuracy of the solver against

When run, the script simulates a time step and evaluates the positions of 4 satellites and estimates the position of the target. A graph is plotted of the target and the satellites in 3D, and diagnostic information is printed to the console.

Kepler Solver
==============
kepler_orbit.m is a script that defines a function kepler_orbit(period,a,b)

The script plots four graphs taken from the solution of the [kepler equation](https://en.wikipedia.org/wiki/Kepler's_equation) along fixed timesteps for one period of orbit

This project was developed in a UC Davis Seminar, Computer Simulations for STEM education
   
