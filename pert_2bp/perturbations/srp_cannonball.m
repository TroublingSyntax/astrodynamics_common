function a_d = srp_cannonball(params)
%SRP_CANNONBALL Computes perturbing SRP acceleration based on 
%   the simplified cannonball model. This function computes the 
%   perturbing acceleration vector in the inertial cartesian frame.
%   Transformation after this function may be required.
%   based on the following parameters and states:
%   srp_params      : parameters for computing srp acceleration, includes:
%   srp_params(1)   : area-to-mass ratio
%   srp_params(2)   : solar flux constant
%   srp_params(3)   : Earth-to-Sun distance
%   srp_params(4:6) : sunlight unit vector

a2mr  = params(1);
G0    = params(2);
d     = params(3);
s_hat = params(4:6);

a_d = a2mr*(G0/d^2)*s_hat;
end

