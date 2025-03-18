function vect_xdot = pert_2bd_keplerian(~,vect_x,args)
%PERT_2BD_KEPLERIAN 1st order perturbed 2BP state function.
%   This function can be passed into ode89 as a function handle
%   to be solved in a keplerian state representation. It expects:
%   time variable t
%   keplerian state vector x
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
a = vect_x(1);
e = vect_x(2);
i = vect_x(3);
RAAN = vect_x(4);
omega = vect_x(5);
p = a*(1-e^2);
h = sqrt(args.mu*p);
b = a*sqrt(1-e^2);

% compute 2 body dynamics
n  = sqrt(args.mu/vect_x(1)^3);
vect_f0 = [0 0 0 0 0 n]';

% Find the equatorial position vector
if vect_x(6) >= 2*pi
    while vect_x(6) >= 2*pi
        vect_x(6) = vect_x(6) - 2*pi;
    end
end
f  = @(f) f - e*sin(f) - vect_x(6);
df = @(f) 1 - e*cos(f);
[~,~,E] = newton_raphson(f,df,vect_x(6),0.0001,10000);
nu = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
rmag = a*(1-e*cos(E));
vect_r_rtn = [rmag 0 0]';
vect_r_eq = RtoE(vect_r_rtn,RAAN,omega+nu,i)';
vect_v_rtn = [] % PICK UP HERE, NEED VEL EQUATORIAL FOR DRAG PARAMS

% update perturbation params with current position info
args.j2_params(4) = rmag;
args.j2_params(5:7) = vect_r_eq;

args.body_params(2:4) = vect_r_eq;

args.drag_params = [args.drag_params; vect_v; rmag];

% find the equitorial perturbing acceleration
vect_a_d_eq = perturbations(...
    drag=args.drag,drag_params=args.drag_params,...
    j2=args.j2,j2_params=args.j2_params,...
    srp=args.srp,srp_params=args.srp_params,...
    body=args.body,body_params=args.body_params);

% convert to rotating frame
vect_a_d_r = EtoR(vect_a_d_eq,RAAN,omega+nu,i);

B=(1/h)*...
[2*a^2*e*sin(nu)                    (2*a^2*p)/rmag              0;
 p*sin(nu)                          (p+rmag)*cos(nu)+rmag*e     0;
 0                                  0                           rmag*cos(nu+omega);
 0                                  0                           (rmag*sin(nu+omega))/sin(i);
 -(p*cos(nu))/e                     ((p+rmag)*sin(nu))/e        -(rmag*sin(nu+omega))/tan(i);
 (b*p*cos(nu))/(a*e)-(2*b*rmag)/a   -(b*(p+rmag)*sin(nu))/(a*e)    0];

% compute xdot including j2 irregular gravity perturbation
vect_xdot = vect_f0 + B*vect_a_d_r;
end

