function [V] = interpmat_nodal2D(x0, y0, p, x1, y1, basis_type)
% [V] = interpmat_nodal2D(x0, y0, x1, y1)
%   compute the interpolation coefficient matrix
%   using Taylor basis
% 
% input: 
%   p:          polynomial order
%   x0, y0:     coordinates of the original data set, length = ((p+1)*(p+2))/2
%   x1, y1:     coordinates of the target data set, of the same length
%   basis_type: type of polynomial basis
% 
% output:
%   V:          interpolation matrix from (x0, y0) to (x1, y1)
%
% WARNING: The floating point stability is important! If one directly uses Taylor basis to evaluate point-values, the error will grow extremely fast with p!!!

% to be done: change the modal basis

N = ((p+1)*(p+2))/2;
m = numel(x1);

[~, ~, ~, ~, Pmn, ~] = matrix_nodal2D(x0, y0, 1.0-x0-y0, p, basis_type);

Tmn = zeros(m, N);  % maps modal to nodal

switch basis_type
case 'Taylor'
    for i = 1: m
        Tmn(i, :) = Taylor_basis_eval(x1(i), y1(i), p);
    end
case 'Bernstein'
    for i = 1: m
        Tmn(i, :) = Bernstein2D_basis_eval(x1(i), y1(i), p);
    end
otherwise
    error('interpmat_nodal2D: unknown polynomial basis type!');
end

V = Tmn / Pmn; % V = Tmn * Pnm;

end