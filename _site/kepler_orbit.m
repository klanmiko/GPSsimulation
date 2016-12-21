function kepler_orbit(period,a,b)
  close all;
  %if the value for the major axis is less than minor axis, swap them
  if a<b
    temp = b;
    b = a;
    a = temp;
  end
  
  e = sqrt(1-b^2/a^2);
  %start at periastron
  xi(1) = a*(1-e);
  yi(1) = 0;
  dt = 0.01;
  
  t = linspace(0,period,period/dt);
  to = 0;
  %plot(xi,yi),
  subplot(2,2,1)
  plot(0,0,'r.');
  %plot(sqrt(a^2-b^2),0,'r.');
  hold on
  %plot(-sqrt(a^2-b^2),0,'r.');
  max_iter = period/dt;
  Eray(1) = 0;
  Mray(1) = 0;
  
  for iterations = 2:max_iter
    M = 2*pi/period*(dt*(iterations-1)-to);
    Mray(iterations) = M;
    E = Eray(iterations-1);
    E = kepler_solve(M,E,e);
    
    xi(iterations) = a*(cos(E)-e);
    yi(iterations) = b*sin(E);
    Eray(iterations) = E;
  end
  
  plot(xi,yi);
  xlabel('semi-major axis (could be a or b)'),ylabel('semi-minor axis'),title("ellipse");
  axis("equal")
  
  subplot(2,2,2);
  plot(t,Eray);
  xlabel('t'),ylabel('E');
  
  subplot(2,2,3);
  plot(t,Mray);
  xlabel('t'),ylabel('M');
  
  subplot(2,2,4);
  plot(Eray,Mray);
  xlabel('E'),ylabel('M');
  printf("eccentricity: %f\n",e);
  
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
