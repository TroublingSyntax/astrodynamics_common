function xdot = two_bd_cartesian(~,x,mu)
%PERT_2BD_CARTESIAN 1st order perturbed 2BP state function.
%   This function can be passed into ode89 as a function handle
%   to be solved in a cartesian state representation. It expects:
%   time variable t
%   cartesian state vector x
%   gravitation parameter of the 2BP mu
%   constant number j2 of spherical harmonics
%   normalizing radius r0

rmag = norm(x(1:3)); % distance
r = x(1:3); % position vector (x,y,z)
v = x(4:6); % velocity vector (vx,vy,vz)

xdot = [v; (-mu/(rmag^3))*r]; % 2 body dynamics
end