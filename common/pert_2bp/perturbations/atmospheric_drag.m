function a_d = atmospheric_drag(params)
%ATMOSPHERIC_DRAG This function computes atmospheric drag near Earth.
%   This function takes in various atmospheric parameters and the
%   spacecraft's velocity in the equitorial frame to compute
%   the acceleration due to atmospheric drag in the equitorial frame.
%   based on the following parameters and states:
%   atmd_params: parameters specific to Earth's atmosphere
%   vect_v_eq: the equitorial velocity vector

% local variables
atmd_params = params(1:6);
vect_v_eq   = params(7:9);
rmag        = params(10);
Cd      = atmd_params(1); % dimensionless
rho0    = atmd_params(2); % kg/km^3
h0      = atmd_params(3); % km
H       = atmd_params(4); % km
r_earth = atmd_params(5); % km
a2mr    = atmd_params(6); % km^2/kg

beta = -(rmag-r_earth-h0)/(H);
v    = norm(vect_v_eq);

a_d = -0.5*Cd*a2mr*rho0*exp(beta)*v^2*(vect_v_eq/v);
end

