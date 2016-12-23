---
# You don't need to edit this file, it's empty on purpose.
# Edit theme's home layout instead if you wanna make some changes
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults
layout: default
---
<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
Synopsis
================
A GPS system allows a receiver to determine its position in xyz coordinates relative to the xyz coordinates of 4 satellites.

The GPS equation:

$$
(x - x_i)^2 + (y - y_i)^2 + (z - z_i)^2 = c^2(t_i - b - s_i)^2
$$

* equates the distance formula with the distance given by travel time
* x<sub>i</sub>, y<sub>i</sub>, z<sub>i</sub> are the coordinates of satellite i
* t<sub>i</sub> is the time message was received in receiver time
* b is the clock bias
* s<sub>i</sub> is the time the message was sent in satellite time

Because there are 4 variables in the function $$f(x,y,z,b)$$ the receiver needs to be within sight of 4 different satellites in order to figure out its position

Method
==========
Using the messages from at least 4 satellites, a GPS solver must then find a solution to a system of nonlinear equations. There exists an algebraic solution for each of the variables, but as per the goal of the course I decided to find a numerical solution for the system of equations. Thus I chose [**Newton's Method**](http://fourier.eng.hmc.edu/e176/lectures/NM/node21.html) to build such a solver.

Simulation
-----------
The issue with building such a solver is that there must be a way to test it. The solver cannot be tested with random values, because the validity and accuracy of the solver cannot be confirmed as the solution will not have been predetermined, and there is no other method I know to solve the GPS equations.

So what is the answer?

The *easier* way is to build a system that can:

1. Take a target coordinate 
2. Calculate the time it takes for a satellite message to travel from a random satellite to the target 
3. Desynchronize the target "clock" and satellite "clock" 
4. Send these parameters to the solver, which solves the system of equations
5. Compare the solution to the target for accuracy.
6. Repeat steps 2-6 over a sample domain

In a real-world scenario, satellites have their own form of GPS that keeps track of their position based on a number of ground stations and an orbital model. To keep the GPS simulator in this project relevant, a decision was made to model orbiting satellites instead of using random data points. The data generated from this system is easier to understand and analyze as compared to random sampling.

However, an orbital model would need to be implemented to simulate the orbit of the satellites.

Kepler Equation
----------------
The [kepler equation](https://en.wikipedia.org/wiki/Kepler's_equation) relates the Mean anomaly (angle) to the eccentric anomaly (another angle) given the eccentricity as a constant (shape of the ellipse)

The eccentric anomaly E is useful because the position of an orbital body can be calculated just knowing this angle.

$$
M = E - esin(E)
$$

This equation has no algebraic solution for E, thus the next best option for solving this equation is again Newton's method. 

The GPS simulator requires a solver for the kepler equation along with the orbital characteristics of each satellite in order to produce their positions at each timestep.

Code Review
================

