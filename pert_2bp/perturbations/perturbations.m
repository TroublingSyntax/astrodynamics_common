function vect_a_d = perturbations(args)
%PERTURBATIONS This function combines many sources of perturbations
%   This function accumulates the acceleration from multiple perturbation
%   sources for ease of use and modularity. takes optional switches
%   for appending different perturbation sources into the return value
%   which is expressed in an equitorial cartesian frame.

arguments
    args.drag = "off";
    args.drag_params = zeros(10,1);
    args.j2   = "off";
    args.j2_params = zeros(7,1);
    args.srp  = "off";
    args.srp_params = zeros(6,1);
    args.body = "off";
    args.body_params = zeros(7,1);
end

vect_a_d = zeros(3,1);

if args.drag == "on"
    vect_a_d = vect_a_d + atmospheric_drag(args.drag_params);
end

if args.j2 == "on"
    vect_a_d = vect_a_d + irreg_grav_j2(args.j2_params);
end

if args.srp == "on"
    vect_a_d = vect_a_d + srp_cannonball(args.srp_params);
end
if args.body == "on"
    vect_a_d = vect_a_d + body(args.body_params);
end

end

