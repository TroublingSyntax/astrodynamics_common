function Ydot = EoMwithDiscreteMatrix(~,Y,uk,nx,A,B,c)
x=Y(1:nx);
Phi=reshape(Y(nx+1:nx+nx^2),nx,nx);
xdot=[x(4);x(5);x(6);uk(1);uk(2);uk(3);-.5086*uk(4)]+c;
Ydot=[xdot;vec(A*Phi);vec(Phi\B);vec(Phi\c)];
end