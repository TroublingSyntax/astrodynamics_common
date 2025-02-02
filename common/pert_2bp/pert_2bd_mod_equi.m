function vect_xdot = pert_2bd_mod_equi(~,vect_x,args)
%PERT_2BD_MOD_EQUI 1st order perturbed 2BP state function.
%   This function can be passed into ode89 as a function handle
%   to be solved in a mod_equi state representation. It expects:
%   time variable t
%   modified equinoctial state vector x
%   gravitation parameter of the 2BP mu

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
p = vect_x(1);
f = vect_x(2);
g = vect_x(3);
h = vect_x(4);
k = vect_x(5);
L = vect_x(6);

RAAN  = atan(k/h);
i     = 2*atan(h/cos(RAAN));
theta = L-RAAN;
omega = atan(g/f)-RAAN;
nu    = theta-omega;
e     = f/cos(omega+RAAN);

rmag = p/(1+e*cos(nu));
vect_r_rtn = [rmag 0 0]';
vect_r_eq  = RtoE(vect_r_rtn,RAAN,theta,i);

% update perturbation params with current position info
args.j2_params(4) = rmag;
args.j2_params(5:7) = vect_r_eq;

args.body_params(2:4) = vect_r_eq;

q = 1+f*cos(L)+g*sin(L);
s_sq = 1+h^2+k^2;


% 2 body dynamics first
vect_f0 = [0 0 0 0 0 sqrt(args.mu*p)*(q/p)^2]';

% perturbations
vect_a_d_eq = perturbations(...
    drag=args.drag,drag_params=args.drag_params,...
    j2=args.j2,j2_params=args.j2_params,...
    srp=args.srp,srp_params=args.srp_params,...
    body=args.body,body_params=args.body_params);

vect_a_d_r  = EtoR(vect_a_d_eq,RAAN,theta,i);

B = sqrt(p/args.mu)*...
    [
     0        (2*p)/q             0;
     sin(L)   ((q+1)*cos(L)+f)/q  -(g*(h*sin(L)-k*cos(L)))/q;
     -cos(L)  ((q+1)*sin(L)+g)/q  (f*(h*sin(L)-k*cos(L)))/q;
     0        0                   (s_sq/(2*q))*cos(L);
     0        0                   (s_sq/(2*q))*sin(L);
     0        0                   (h*sin(L)-k*cos(L))/q
     ];

vect_xdot = vect_f0 + B*vect_a_d_r;
end

