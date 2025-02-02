function vect_xdot = pert_2bd_milankovitch(~,vect_x,args)
%PERT_2BD_MOD_EQUI 1st order perturbed 2BP state function.
%   This function can be passed into ode89 as a function handle
%   to be solved in a mod_equi state representation. It expects:
%   time variable t
%   modified equinoctial state vector x
%   gravitation parameter of the 2BP mu
%   constant number j2 of spherical harmonics
%   normalizing radius r0
%   srp_params parameters for solar radiation pressure perturbation

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

% local variables
vect_h = vect_x(1:3);
vect_e = vect_x(4:6);
L      = vect_x(7);

h = norm(vect_h);
e = norm(vect_e);
p = h^2/args.mu;
a = p/(1-e^2);

vect_h_hat = vect_h/h;
vect_e_hat = vect_e/e;

% compute n_hat in inertial frame
hx = vect_h_hat(1);
hy = vect_h_hat(2);
ny = sqrt(hx^2/(hy^2+hx^2));
nx = -(ny*hy)/hx;
vect_n_hat = [nx ny 0]';
if dot(vect_h_hat,vect_n_hat) ~= 0
    % handle the sqrt taking away the negative option for ny
    ny = -sqrt(hx^2/(hy^2+hx^2));
    nx = -(ny*hy)/hx;
    vect_n_hat = [nx ny 0]';
    if dot(vect_h_hat,vect_n_hat) ~= 0
        % if problem isn't fixed, there is problem
        warning("computing node vector failed!")
    end
end

RAAN  = atan(ny/nx);
omega = acos(dot(vect_n_hat,vect_e_hat));
nu    = L - omega - RAAN;
i = acos(vect_h_hat(3));

rmag = p/(1+e*cos(nu));
vect_r_r = [rmag 0 0]';
vect_r_eq = RtoE(vect_r_r,RAAN,nu+omega,i);

% update perturbation params with current position info
args.j2_params(4) = rmag;
args.j2_params(5:7) = vect_r_eq;

v = sqrt(args.mu*((2/rmag)-(1/a)));
fpa = acos(sqrt(args.mu*p)/(rmag*v));

% quadrant check for flight path angle
if nu > pi
    if fpa > 0
        fpa = -fpa;
    end
elseif nu < pi
    if fpa < 0
        fpa = -fpa;
    end
end

vect_v_r = [v*sin(fpa) v*cos(fpa) 0]';
vect_h_r = EtoR(vect_h,RAAN,nu+omega,i);

% 2 body dynamics first
vect_f0 = [0 0 0 0 0 0 h/rmag^2]';

% perturbations
vect_a_d_eq = perturbations(...
    drag=args.drag,drag_params=args.drag_params,...
    j2=args.j2,j2_params=args.j2_params,...
    srp=args.srp,srp_params=args.srp_params,...
    body=args.body,body_params=args.body_params);

B = [skewSymmetric(RtoE(vect_r_r,RAAN,nu+omega,i));
     (1/args.mu)*(skewSymmetric(RtoE(vect_v_r,RAAN,nu+omega,i))*skewSymmetric(RtoE(vect_r_r,RAAN,nu+omega,i))-skewSymmetric(vect_h));
     ((vect_r_eq(3))/(h*(h+vect_h(3))))*RtoE(vect_h_r,RAAN,nu+omega,i)'];

vect_xdot = vect_f0 + B*vect_a_d_eq;
end