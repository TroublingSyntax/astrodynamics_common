function xdot = two_bd_cartesian_model(~,x,xi)
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

xdot = [xi(1)*v(1) xi(2)*v(2) xi(3)*v(3) ...
    xi(4)*r(1)*rmag^-3 xi(5)*r(2)*rmag^-3 ...
    xi(6)*r(3)*rmag^-3]'; % 2 body dynamics
end