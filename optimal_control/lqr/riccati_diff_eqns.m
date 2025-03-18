function Kdot = riccati_diff_eqns(~,K,A,B,Q,R,P)
%RICCATI_DIFF_EQNS Riccati differential equations
%   This function is to be numerically integrated to solve the optimal gain
%   matrix K for an LQR controller based on the following input arguments:
%   A - system matrix
%   B - control input matrix
%   Q - weight matrix for state error
%   R - weight matrix for control input
%   P - weight matrix for both state and control input
%   It returns the current time derivative of the gain matrix K:
%   Kdot - time derivative of gain matrix K
arguments
    ~
    K
    A
    B
    Q
    R
    P
end
%addpath(genpath('../../common/utils'))
K=reshape(K,[4 4]);
Kdot=vec((K*B/R+P/R)*(P'+B'*K)-K*A-Q-A'*K);
end