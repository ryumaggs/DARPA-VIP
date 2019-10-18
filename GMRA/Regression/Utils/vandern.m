function A = vandern(v,n)

%VANDER Vandermonde matrix.
%   A = VANDER(V) returns the Vandermonde matrix whose columns
%   are powers of the vector V, that is A(i,j) = v(i)^(n-j).
%
%   Class support for input V:
%      float: double, single

%   Small modification of Matlab's vander

v = v(:);
m = length(v);
if m == 0
    A = reshape(v, m, m);
    return
end
A = repmat(v, 1, n);
A(:, 1) = 1;
A = cumprod(A, 2);

return