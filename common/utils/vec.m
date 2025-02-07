function x = vec(M)
%VEC Vectorizes a matrix
%   This function takes a matrix as an input and sequentially converts it's
%   columns to a column vector.
m = size(M,1);
n = size(M,2);
x = zeros(m*n,1);
j = 0;
for i=1:n
    x(j*m+1:j*m+m,1) = M(:,i);
    j = j+1;
end
end

