function vect_xdot = pert_2bd_cartesian(~,vect_x,args)
%PERT_2BD_CARTESIAN 1st order perturbed 2BP state function.
%   This function can be passed into ode89 as a function handle
%   to be solved in a cartesian state representation. It expects:
%   time variable t
%   cartesian state vector x

arguments
    ~
    vect_x
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


addpath(genpath('perturbations'))

rmag = norm(vect_x(1:3)); % distance
vect_r = vect_x(1:3); % position vector (x,y,z)
vect_v = vect_x(4:6); % velocity vector (vx,vy,vz)

% update perturbation params with current position info
args.j2_params(4) = rmag;
args.j2_params(5:7) = vect_r;

args.body_params(2:4) = vect_r;

args.drag_params = [args.drag_params; vect_v; rmag];

vect_f0  = [vect_v; (-args.mu/(rmag^3))*vect_r]; % 2 body dynamics
B   = [zeros(3); eye(3)]; % mapping of perturbation into state

vect_a_d = perturbations(...
    drag=args.drag,drag_params=args.drag_params,...
    j2=args.j2,j2_params=args.j2_params,...
    srp=args.srp,srp_params=args.srp_params,...
    body=args.body,body_params=args.body_params);

vect_xdot = vect_f0 + B*vect_a_d;
end