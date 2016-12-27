---
# You don't need to edit this file, it's empty on purpose.
# Edit theme's home layout instead if you wanna make some changes
# See: https://jekyllrb.com/docs/themes/#overriding-theme-defaults
layout: default
---
<script src="https://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>
Demo
===========
<blockquote class="imgur-embed-pub" lang="en" data-id="75mQz3V"><a href="//imgur.com/75mQz3V"></a></blockquote><script async src="//s.imgur.com/min/embed.js" charset="utf-8"></script>

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

Kepler Orbit
------------
The kepler solver is described in the file kepler_orbit.m. The functional portions of the code will be explained as follows.

~~~~
function kepler_orbit(period,a,b)
...
e = sqrt(1-b^2/a^2);
~~~~

The eccentricity of the elliptical orbit is determined from the major and minor axis.

~~~~
dt = 0.01;
t = linspace(0,period,period/dt);
max_iter = period/dt;
Eray(1)=0;
Mray(1)=0;
~~~~

t is a vector from 0 to the period spaced evenly by dt, Eray and Mray hold the values of the eccentric anomaly and mean anomaly at each time step.

~~~~
for iterations = 2:max_iter
  M = 2*pi/period*(dt*(iterations-1)-to);
    Mray(iterations) = M;
    E = Eray(iterations-1);
    E = kepler_solve(M,E,e);
    
    xi(iterations) = a*(cos(E)-e);
    yi(iterations) = b*sin(E);
    Eray(iterations) = E;
end
~~~~

The script iterates over the entire ellipse, calculating the value of E and the coordinates of the orbital body at each timestep. The last value of E is used as the initial guess for the current value.

~~~~
function out = kepler(E,e)
         out = E - e*sin(E);
         
function out = pkepler(E,e)
         out = 1 - e*cos(E);
         
function out = kepler_solve(M,Eo,e) %Eo is for initial guess
  count = 0;
  internal_max = 20000;
  E = Eo;
  
  while abs(kepler(E,e)-M) > 10e-8 && count<internal_max
    E = E - (kepler(E,e)-M)/pkepler(E,e);
    count++;
  end  
  out = E;
~~~~
The three function above define the kepler equation, its derivative, and an application of newton's method to solve the equation.

For an eccentricity of 0.661438, the following plot is obtained.
![kepler]({{site.url}}/images/r7f1kQS.png)

GPS Solver
--------------
The code for the GPS solver is complex and contains parts of the above kepler solver in its implementation, which is explained below.

~~~~
function GPS(target) % target = [0,-1,0] for example
disp(target);
close all;
%note z-value cannot be zero or the jacobian is singular
%periastron time is the time when the satellite passes the periastron
satellite1 = [30,4,4,60,3]; %[period a b periastron_time z]
satellite2 = [100,4,4,30,0.3];
satellite3 = [100,4,4,0,0.7];
satellite4 = [50,4,4,7,-1];
c = 10;
bias = 20;
[XS,YS,ZS] = sphere(50);
current_time = 0;
dt = 1;
~~~~

The initial conditions for the satellites are set up, the speed of light is set to 10, the desynchronization of satellite and target clocks is set to 20, and the time step is set to 1.

~~~~
while current_time < 200
  
  ps1 = kepler_orbit(satellite1(1),satellite1(2),satellite1(3),satellite1(4),current_time);
  ps2 = kepler_orbit(satellite2(1),satellite2(2),satellite2(3),satellite2(4),current_time);
  ps3 = kepler_orbit(satellite3(1),satellite3(2),satellite3(3),satellite3(4),current_time);
  ps4 = kepler_orbit(satellite4(1),satellite4(2),satellite4(3),satellite4(4),current_time);
~~~~

Here ps<sub>i</sub> is a vector representing the $$(x,y)$$ coordinates of each satellite determined from its orbit using the kepler solver.

~~~~  
  t = zeros(1,4);
  
  t(1) = findTrueReceiveTimeForSatellite(target,[ps1(1),ps1(2),satellite1(5)],c)+current_time; % target, [x,y,z], c
  t(2) = findTrueReceiveTimeForSatellite(target,[ps2(1),ps2(2),satellite2(5)],c)+current_time;
  t(3) = findTrueReceiveTimeForSatellite(target,[ps3(1),ps3(2),satellite3(5)],c)+current_time;
  t(4) = findTrueReceiveTimeForSatellite(target,[ps4(1),ps4(2),satellite4(5)],c)+current_time;
~~~~

The t vector holds the value of the simulated travel times between each satellite and the target, adjusted to the timestep.

~~~~  
  x = trilaterate(c,[ps1(1),ps1(2),satellite1(5),t(1)+bias],
  [ps2(1),ps2(2),satellite2(5),t(2)+bias],
  [ps3(1),ps3(2),satellite3(5),t(3)+bias],
  [ps4(1),ps4(2),satellite4(5),t(4)+bias],
  current_time);
~~~~

Here x is a vector representing the coordinates of the target given by the GPS solver. The clock desynchronization is added to the t vector before it is passed to the solver.

~~~~
function out = findTrueReceiveTimeForSatellite(target,psatellite,c)
  out = sqrt((target(1)-psatellite(1))^2+(target(2)-psatellite(2))^2+(target(3)-psatellite(3))^2)/c;
end
~~~~

The 	findTrueReceiveTimeForSatellite	function simply uses the distance formula to find the time light takes to travel from each satellite to the target.

~~~~
function out = trilaterate(c,d1,d2,d3,d4,sent_time)
  x = zeros(4,1);
  c
  sent_time
  d1
  d2
  d3
  d4
  x(1) = 10; x(2) = 10; x(3) = 10; x(4) = 3; %x = [x y z t]

  TOL = 1e-10; error = 1;
  count = 0;
~~~~

The function 	trillaterate	takes as arguments the speed of light, coordinates, message receive time in target frame, and the send time for all of the messages (assuming the satellites all send messages at the same time). It defines the x vector to hold the x,y,z coordinates of the target to be iterated upon, and x(4) represents the clock bias. 

~~~~
while max(abs(error)) > TOL
xold = x;

error = [(x(1)-d1(1))^2+(x(2)-d1(2))^2+(x(3)-d1(3))^2-c^2*(-x(4)+d1(4)-sent_time)^2;
(x(1)-d2(1))^2+(x(2)-d2(2))^2+(x(3)-d2(3))^2-c^2*(-x(4)+d2(4)-sent_time)^2;
(x(1)-d3(1))^2+(x(2)-d3(2))^2+(x(3)-d3(3))^2-c^2*(-x(4)+d3(4)-sent_time)^2;
(x(1)-d4(1))^2+(x(2)-d4(2))^2+(x(3)-d4(3))^2-c^2*(-x(4)+d4(4)-sent_time)^2];

jacobian = [2*(x(1)-d1(1)), 2*(x(2)-d1(2)), 2*(x(3)-d1(3)), +c^2*2*(-x(4)+d1(4)-sent_time);
2*(x(1)-d2(1)), 2*(x(2)-d2(2)), 2*(x(3)-d2(3)), +c^2*2*(-x(4)+d2(4)-sent_time);
2*(x(1)-d3(1)), 2*(x(2)-d3(2)), 2*(x(3)-d3(3)), +c^2*2*(-x(4)+d3(4)-sent_time);
2*(x(1)-d4(1)), 2*(x(2)-d4(2)), 2*(x(3)-d4(3)), +c^2*2*(-x(4)+d4(4)-sent_time)]

if(det(jacobian)!=0) x = x-inv(jacobian)*error;
else 
  display("jacobian is singular error");
  break;
endif
error
[count x(1) x(2) x(3) x(4)]
count+=1;
endwhile
   r = [c^2*(-x(4)+d1(4)-sent_time)^2;
 c^2*(-x(4)+d2(4)-sent_time)^2;
 c^2*(-x(4)+d3(4)-sent_time)^2;
 c^2*(-x(4)+d4(4)-sent_time)^2].^(1/2)
 
   z = [-x(4)+d1(4);
  -x(4)+d2(4);
  -x(4)+d3(4);
  -x(4)+d4(4)]
  
  out = x;
end
~~~~

The error function is the GPS equation with the right hand side subtracted to give an equation that equals zero. Then the jacobian is determined using the partial derivatives of the GPS equation with respect to x,y,z and bias, each row representing the data from each satellite. 

The function checks that the jacobian is not singular, and then performs newton's method to approximate the coordinates of the target. The r vector gives the distance from each satellite to the target in terms of time for each message to be received. The z vector gives the true receive time for each message (in satellite frame).
	
Conclusion
=============
From running the simulation to completion a number of times, the simulator has shown to be effective most of the time on a qualitative basis. There are cases in which the solver converges to the wrong root for certain positions of the satellites. The error rate also seems to increase with the scaling of the speed of light, where values in the extremes cause the solver to error wildly. It is beyond the scope of this project to test or find solutions to the errors in the GPS solver. However, the easiest solution would be to use previously found values in the solver as the initial guess for the soler(assuming that the target does not move very far). An in-depth study may focus on producing a fractal of which initial guesses produce the correct solution for the GPS solver. Another system would have to be developed to correctly scale satellite coordinates and the value of the speed of light, in a real-world situation the time travelled between satellite and receiver is significantly smaller than the distance between the satellite and receiver, thus a way to scale these two values would need to be found while maintaining the accuracy of the GPS system.

Another limitation is that the kepler solver only produces x and y coordinates, when in reality satellites orbit at an angle, and the current solution does not account for this tilt. This issue is further affected by the trillateration function's inability to deal with a singular matrix, thus forcing the implementation to create fake z "levels" for each orbital ellipse. 

Thanks to Dr. Hafez and Dr. Tavernetti of UC Davis for helping me with this project.
