function [M, Minv, Dx, Dy, Pmn, Pnm] = matrix_nodal2D(x, y, z, p)
% [M, Minv, Dx, Dy, Pmn, Pnm] = matrix_nodal2D(x, y, z, p)
%   compute the matrices for p-th order polynomial space on the right triangle
%   using Taylor basis
% 
% input: 
%   x,y,z:  barycentric coordinates of collocation nodes, length = ((p+1)*(p+2))/2
%   p:      polynomial order >= 0
% 
% output: 
%   M:      mass matrix with respect to nodal DOFs, size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]
%   Minv:   inverse of the mass matrix, size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]
%   Dx, Dy: differentiation matrix, size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]
%   Pmn:    maps modal coefs to nodal values (accurate), size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]
%   Pnm:    maps nodal values to modal coefs (should be avoided as possible), size = [((p+1)*(p+2))/2, ((p+1)*(p+2))/2]

% checked.

id2D = @(i,j,p) ((2*p+4-i)*(i-1))/2 + j;
N = ((p+1)*(p+2))/2;

% Currently, we use the Taylor basis for construction.

% Pmn maps modal coefs to nodal values (accurate)
Pmn = zeros(N, N);
for i = 1: p+1
    for j = 1: p+2-i
        Pmn(id2D(i,j,p), :) = Taylor_basis_eval(x(id2D(i,j,p)), y(id2D(i,j,p)), p);
    end
end

% Pnm maps nodal values to modal coefs (big floating-point error, should be avoided as possible)
Pnm = inv(Pmn);

% compute modal matrices
[M_m, Dx_m, Dy_m] = Taylor_basis_mat(p);

% compute nodal matrices
Dx = Pmn * Dx_m / Pmn; % Dx = Pmn * Dx_m * Pnm;
Dy = Pmn * Dy_m / Pmn; % Dy = Pmn * Dy_m * Pnm;

% M = Pnm' * M_m * Pnm;
M = Pmn' \ M_m / Pmn; % it proves to be more accurate by numerical experiments
M = 0.5*(M + M');

% Minv = Pmn * inv(M_m) * Pmn';
Minv = inv(M); % it proves to be more accurate by numerical experiments
Minv = 0.5*(Minv + Minv');


end