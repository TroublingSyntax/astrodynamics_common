function [dydt] = EoM(t,y,A,B,Q,R,P,ktraj)
%EOM Equations of Motion of LQR system
%   This function is to be numerically integrated to solve for the
%   augmented state equation y. The inputs include
%   t the time
%   y the augmented state
%   A - system matrix
%   B - control input matrix
%   Q - weight matrix for state error
%   R - weight matrix for control input
%   P - weight matrix for both state and control input
%   ktraj - optimal gain matrix over time expressed as a vector
arguments
    t
    y
    A
    B
    Q
    R
    P
    ktraj
end
x=y(1:4,1);
lambda=y(5:8,1);
K=deval(ktraj,t);
K=reshape(K,[4 4]);
S=-inv(R)*P'-R\B'*K;
u=S*x;
xdot=(A+B*S)*x;
lambdadot=-Q*x-P*u-A'*lambda;
dydt=[xdot;lambdadot;u];
end

