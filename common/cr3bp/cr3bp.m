function vect_xdot = cr3bp(~,vect_x,mu)
%CR3BP ode propagation function for circular restricted 3 body problem.
%   The cr3bp function takes non-dimensional cartesian state information
%   x and the mass ratio mu as inputs. It returns xdot using equations of
%   motion also in non-dimensional units. The state information is in the
%   synodic frame of the two significant masses of the problem.

% local variables
x    = vect_x(1);
y    = vect_x(2);
z    = vect_x(3);
xdot = vect_x(4);
ydot = vect_x(5);
zdot = vect_x(6);

r13 = sqrt((x+mu)^2 + y^2 + z^2);
r23 = sqrt((x-1+mu)^2 + y^2 + z^2);

% equations of motion
xddot = 2*ydot + x - (((1-mu)*(x+mu))/r13^3) - ((mu*(x-1+mu))/r23^3);
yddot = -2*xdot + y - (((1-mu)*y)/r13^3) - ((mu*y)/r23^3);
zddot = -(((1-mu)*z)/r13^3) - ((mu*z)/r23^3);

% return the xdot vector
vect_xdot = [xdot ydot zdot xddot yddot zddot]';

end