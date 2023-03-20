function [V] = interpmat_nodal2D(x0, y0, p, x1, y1)
% [V] = interpmat_nodal2D(x0, y0, x1, y1)
%   compute the interpolation coefficient matrix
%   using Taylor basis
% 
% input: 
%   p:          polynomial order
%   x0, y0:     coordinates of the original data set, length = ((p+1)*(p+2))/2
%   x1, y1:     coordinates of the target data set, of the same length
% 
% output:
%   V:          interpolation matrix from (x0, y0) to (x1, y1)

N = ((p+1)*(p+2))/2;
m = numel(x1);

[~, ~, ~, ~, Pnm] = matrix_nodal2D(x0, y0, 1.0-x0-y0, p);
Tmn = zeros(m, N);  % maps modal to nodal
for i = 1: m
    Tmn(i, :) = Taylor_basis_eval(x1(i), y1(i), p);
end
V = Tmn * Pnm;

end