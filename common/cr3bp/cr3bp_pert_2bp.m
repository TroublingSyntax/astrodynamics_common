function xdot = cr3bp_pert_2bp(~,x,args)
%CR3BP_PERT_2BP Circular restricted 3 body problem expressed in 2bp.
%   This function works very similarly to the normal perturbed 2bp
%   function, but this also incorporates extra state information so that
%   the position of the second smaller mass can also have it's position
%   propagated over time. it assumes that the smaller mass has a circular
%   orbit around the larger, and that mass 3 experiences the gravity of
%   mass 2 as a perturbation in a 2bp with mass 1.

arguments
    ~
    x
    args.mu_mass1
    args.mu_m1m2
    args.mu
    args.drag = "off"
    args.drag_params = zeros(10,1);
    args.j2 = "off"
    args.j2_params = zeros(7,1);
    args.srp = "off"
    args.srp_params = zeros(6,1);
    args.body = "off";
    args.body_params = zeros(7,1);
end


addpath(genpath('../pert_2bp/perturbations'))

rmag = norm(x(1:3)); % distance
r = x(1:3); % position vector (x,y,z)
v = x(4:6); % velocity vector (vx,vy,vz)

r_m2_mag = norm(x(7:9)); % m2 distance
r_m2     = x(7:9); % position vector (x,y,z)
v_m2     = x(10:12); % velocity vector (vx,vy,vz)

f0  = [v; (-args.mu_mass1/(rmag^3))*r; ...
       v_m2; (-args.mu_m1m2/(r_m2_mag^3))*r_m2]; % 2 body dynamics
B   = [zeros(3); eye(3); zeros(6,3)]; % mapping of perturbation into state

args.body_params(2:4) = r;
args.body_params(5:7) = r_m2;

args.j2_params(4) = rmag;
args.j2_params(5:7) = r;

a_d = perturbations(drag=args.drag,drag_params=args.drag_params,...
    j2=args.j2,j2_params=args.j2_params,srp=args.srp,...
    srp_params=args.srp_params,body=args.body,body_params=args.body_params);

xdot = f0 + B*a_d;
end

