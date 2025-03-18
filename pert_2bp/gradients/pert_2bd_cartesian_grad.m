function dfdx = pert_2bd_cartesian_grad(t,x,args)
arguments
    t
    x
    args.j2 = "off"
    args.j2_params = 0;
    args.srp = "off"
    args.srp_params = 0;
end
% local variables
mu = args.j2_params(1);
j2 = args.j2_params(2);
r0 = args.j2_params(3);
a2mr  = args.srp_params(1);
G0    = args.srp_params(2);
d     = args.srp_params(3);

r = x(1:3,1);
x = r(1);
y = r(2);
z = r(3);
rmag = norm(r);
rhatt = (r/rmag)';
n1t = [1 0 0];
n2t = [0 1 0];
n3t = [0 0 1];

da_j2 = (-mu/rmag^3)*eye(3)+((3*mu)/rmag^4)*r*rhatt-((3*mu*j2*r0^2)/(2*rmag^5))*...
[((-5*x)/(rmag^2))*(1-(7*z^2)/(rmag^2))*r'+(1-(5*z^2)/(rmag^2))*n1t+((10*x*z)/(rmag^2))*n3t;
 ((-5*y)/(rmag^2))*(1-(7*z^2)/(rmag^2))*r'+(1-(5*z^2)/(rmag^2))*n2t-((10*y*z)/(rmag^2))*n3t;
 ((5*z)/(rmag^2))*(((7*z^2)/(rmag^2))-3)*r'+3*(1-(5*z^2)/(rmag^2))*n3t];

n = 0.00088233532457904;

s = r - d*([cos(n*t) sin(n*t) 0]');
smag = norm(s);

da_srp = a2mr*(G0/d^2)*(-smag^-3*(s*s')+(smag^-1)*eye(3));

da = da_j2 + da_srp;

dfdx = ...
    [zeros(3) eye(3);
     da       zeros(3)];
end