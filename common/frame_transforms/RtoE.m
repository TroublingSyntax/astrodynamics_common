function v_e = RtoE(v_r,omega,theta,i)
%RtoP This function transforms a vector from the
%     rotating frame (r,theta,h) to the equitorial frame.
%   INPUT: the vector given in rotating frame components,
%          the right ascension of the ascending node, theta
%          and the inclination. these angles need to be in radians
%
%   OUTPUT: the original vector in terms of the equitorial frame.

% a11
a=(cos(omega)*cos(theta))-(sin(omega)*cos(i)*sin(theta));
% a12
b=(-cos(omega)*sin(theta))-(sin(omega)*cos(i)*cos(theta));
% a13
c=sin(omega)*sin(i);
% a21
d=(sin(omega)*cos(theta))+(cos(omega)*cos(i)*sin(theta));
% a22
e=(-sin(omega)*sin(theta))+(cos(omega)*cos(i)*cos(theta));
% a23
f=-cos(omega)*sin(i);
% a31
g=sin(i)*sin(theta);
% a32
h=sin(i)*cos(theta);
% a33
i=cos(i);

F=[a b c;
   d e f;
   g h i];

v_e = F * v_r;
end