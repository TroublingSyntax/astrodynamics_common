function [vect_t_dim,vect_x_dim] = dimensionalize(X,vect_t,l,t)
%DIMENSIONALIZE Dimensionalizes non-dimensional cartesian cr3bp state
%   This function takes in a cartesian, non-dimensional state vector
%   information as a time series matrix, and dimensionalizes specified
%   by the following input parameters:
%   x: cartesian state information (non-dimensional)
%   l: characteristic length
%   t: characteristic time

% local variables
vect_r_rtn = X(:,1:3);
vect_v_rtn = X(:,4:6);

% dimensionalize position vector
vect_r_rtn_dim = l * vect_r_rtn;

% dimensionalize time
vect_t_dim = t * vect_t;

% dimensionalize velocty vector
vect_v_rtn_dim = (l/t) * vect_v_rtn;

vect_x_dim = [vect_r_rtn_dim vect_v_rtn_dim];
end

