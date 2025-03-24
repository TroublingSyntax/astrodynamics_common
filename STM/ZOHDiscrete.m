function [Ak,Bk,ck,Yk1] = ZOHDiscrete(tk,tk1,Yk,A,B,c,nx,nu,opts)
traj=ode45(@(t,Y)...
EoMwithDiscreteMatrix(t,Y,zeros(nu,1),nx,A,B,c),[tk tk1],Yk,opts);
Yk1=deval(traj,tk1);
Ak=reshape(Yk1(nx+1:nx+nx^2),nx,nx);
Bk=Ak*[reshape(Yk1(nx+nx^2+1:nx+nx^2+nx*nu),nx,nu)];
ck=Ak*[Yk1(nx+nx^2+nx*nu+1:nx+nx^2+nx*nu+nx)];
end

