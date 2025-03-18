function rtilde = skewSymmetric(r)
%SKEWSYMMETRIC Creates skew-symmetric matrix from 3d vector
%   This function takes a single 3x1 vector as an input
%   and produces a skew-symmetric matrix that, when pre-multiplied
%   to another vector, performs the same matrix operation as a vector
%   cross product. For example, if r is the input vector, then rtilde
%   is the output of this function, and rtilde*v = cross(r,v) = the
%   cross product of r and v (r x v).
rtilde = [0 -r(3) r(2);
          r(3) 0 -r(1);
          -r(2) r(1) 0];
end

