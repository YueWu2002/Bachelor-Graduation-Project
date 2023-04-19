function [M, Minv, D, Pmn, Pnm, Tl, Tr] = matrix_nodal1D(x)
% [M, Minv, D, Pmn, Pnm, Tl, Tr] = matrix_nodal1D(x)
%   compute the matrices for p-th order polynomial space on [-1,1]
%   using Legendre basis
% 
% input: 
%   x:  coordinate of collocation nodes, length = p+1
% 
% output: 
%   M:      mass matrix with respect to nodal DOFs, size = [p+1, p+1]
%   Minv:   inverse of the mass matrix, size = [p+1, p+1] (exact)
%   D:      differentiation matrix, size = [p+1, p+1]
%   Pmn:    maps modal coefs to nodal values, size = [p+1, p+1]
%   Pnm:    maps nodal values to modal coefs, size = [p+1, p+1]
%   Tl:     maps nodal values to left trace, size = [1, p+1]
%   Tr:     maps nodal values to right trace, size = [1, p+1]

% checked.

p = length(x) - 1;

% Pmn maps modal coefs to nodal values (exact)
Pmn = zeros(p+1, p+1);
for i = 1: p+1
    Pmn(:, i) = legendrep(x, i-1);
end

% Pnm maps nodal values to modal coefs, size = [p+1, p+1]
Pnm = inv(Pmn);

% compute modal matrices (exact)
[M_m, Minv_m, D_m, Tl_m, Tr_m] = Legendre_basis_mat(p);

% compute nodal matrices
D = Pmn * D_m / Pmn; % D = Pmn * D_m * Pnm;

Minv = Pmn * Minv_m * Pmn';
Minv = 0.5*(Minv + Minv');

M = inv(Minv);  % here, we use inv() since Minv is exact while M is not
% M = Pnm' * M_m * Pnm;
M = 0.5*(M + M');

Tl = interpolatemat(x, -1);
Tr = interpolatemat(x, 1);

end