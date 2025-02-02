function vect_a_d = body(params)
%BODY Computes the perturbing acceleration of an additional body.
%   This function computes the perturbing acceleration on a standard
%   two body problem. It takes as inputs:
%   mu_body    : the gravitational parameters of the extra body
%   vect_r     : the cartesian position vector of the object of interest
%   vect_r_body: the cartesian position vector of the perturbing body

mu_body = params(1);
vect_r = params(2:4);
vect_r_body = params(5:7);

vect_a_d = -mu_body*(((vect_r-vect_r_body)/(norm(vect_r-vect_r_body)^3))...
         + ((vect_r_body)/(norm(vect_r_body)^3)));
end

