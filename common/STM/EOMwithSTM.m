function Xdot = EOMwithSTM(t,X,nx,f,dfdx,args)
arguments
    t
    X
    nx
    f
    dfdx
    args.j2 = "off";
    args.j2_params = 0;
    args.srp = "off";
    args.srp_params = 0;
end
x = X(1:nx);
phi = reshape(X(nx+1:nx+nx^2),nx,nx);
A = dfdx(t,x,j2=args.j2,j2_params=args.j2_params,...
    srp=args.srp,srp_params=args.srp_params);
Xdot = [f(t,x,mu=args.j2_params(1),j2=args.j2,j2_params=args.j2_params,...
    srp=args.srp,srp_params=args.srp_params);
        vec(A*phi)];
end

