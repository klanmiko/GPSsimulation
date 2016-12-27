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
while current_time < 200
  
  ps1 = kepler_orbit(satellite1(1),satellite1(2),satellite1(3),satellite1(4),current_time);
  ps2 = kepler_orbit(satellite2(1),satellite2(2),satellite2(3),satellite2(4),current_time);
  ps3 = kepler_orbit(satellite3(1),satellite3(2),satellite3(3),satellite3(4),current_time);
  ps4 = kepler_orbit(satellite4(1),satellite4(2),satellite4(3),satellite4(4),current_time);
  
  t = zeros(1,4);
  
  t(1) = findTrueReceiveTimeForSatellite(target,[ps1(1),ps1(2),satellite1(5)],c)+current_time; % target, [x,y,z], c
  t(2) = findTrueReceiveTimeForSatellite(target,[ps2(1),ps2(2),satellite2(5)],c)+current_time;
  t(3) = findTrueReceiveTimeForSatellite(target,[ps3(1),ps3(2),satellite3(5)],c)+current_time;
  t(4) = findTrueReceiveTimeForSatellite(target,[ps4(1),ps4(2),satellite4(5)],c)+current_time;
  
  x = trilaterate(c,[ps1(1),ps1(2),satellite1(5),t(1)+bias],
  [ps2(1),ps2(2),satellite2(5),t(2)+bias],
  [ps3(1),ps3(2),satellite3(5),t(3)+bias],
  [ps4(1),ps4(2),satellite4(5),t(4)+bias],
  current_time);
  
  subplot(1,2,1);
  surf(XS,YS,ZS), shading flat;
  hold on
  plot3(ps1(1),ps1(2),satellite1(5),'.','markersize',20),
  plot3(ps2(1),ps2(2),satellite2(5),'.','markersize',20),
  plot3(ps3(1),ps3(2),satellite3(5),'.','markersize',20),
  plot3(ps4(1),ps4(2),satellite4(5),'.','markersize',20);
  plot3(x(1),x(2),x(3),'r.','markersize',20);
  axis([-10,10,-10,10,-10,10]);
  view(2);
  hold off
  
  subplot(1,2,2);
  surf(XS,YS,ZS), shading flat;
  hold on
  plot3(ps1(1),ps1(2),satellite1(5),'.','markersize',20),
  plot3(ps2(1),ps2(2),satellite2(5),'.','markersize',20),
  plot3(ps3(1),ps3(2),satellite3(5),'.','markersize',20),
  plot3(ps4(1),ps4(2),satellite4(5),'.','markersize',20);
  plot3(x(1),x(2),x(3),'r.','markersize',20);
  axis([-10,10,-10,10,-10,10]);
  view(3);
  hold off
  
  disp("should be zero");
  disp(x(4)-bias);
  disp(x(1)-target(1));
  disp(x(2)-target(2));
  disp(x(3)-target(3));
  pause(1);
  current_time += dt;
end
end
function out = findTrueReceiveTimeForSatellite(target,psatellite,c)
  out = sqrt((target(1)-psatellite(1))^2+(target(2)-psatellite(2))^2+(target(3)-psatellite(3))^2)/c;
end
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
function out = kepler_orbit(period,a,b,periastron_time,current_time)
  if a<b
    temp = b;
    b = a;
    a = temp;
  end
  
  e = sqrt(1-b^2/a^2);
  M = 2*pi/period*(current_time-periastron_time);
  E = kepler_solve(M,pi,e);
  
  out = [a*(cos(E)-e),b*sin(E)];
end
function out = kepler(E,e)
         out = E - e*sin(E);
end
function out = pkepler(E,e)
         out = 1 - e*cos(E);
end
function out = kepler_solve(M,Eo,e) %Eo is for initial guess
  count = 0;
  internal_max = 20000;
  E = Eo;
  while abs(kepler(E,e)-M) > 10e-8 && count<internal_max
    E = E - (kepler(E,e)-M)/pkepler(E,e);
    count++;
  end
  out = E;
end
