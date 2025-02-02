function a_d = irreg_grav_j2(params)
%IRREG_GRAV_J2 Computes perturbing j2 irregular gravity acceleration.
%   This function computes the perturbing acceleration vector in inertial
%   cartesian frame. Transformation after this function may be required.
%   based on the following parameters and states:
%   mu  : the gravitational parameter
%   j2  : body j2 value
%   r0  : normalizing radius
%   rmag: magnitude of position vector
%   r   : position vector

% local variables
mu = params(1);
j2 = params(2);
r0 = params(3);
rmag = params(4);
r = params(5:7);

factor = (-3*mu*j2*r0^2)/(2*rmag^5);
a_d = factor*[(1-5*(r(3)^2/rmag^2))*r(1); ...
    (1-5*(r(3)^2/rmag^2))*r(2); ...
    (3-5*(r(3)^2/rmag^2))*r(3)];
end