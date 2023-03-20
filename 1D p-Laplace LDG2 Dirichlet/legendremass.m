function [M, Minv] = legendremass(n)
% [M, Minv] = legendremass(n)
%   compute the mass matrix and its inverse for Legendre polynomials of order 0---n in [-1,1]
%
% input: 
%   n:  max order of Legendre polynomials
% 
% output: 
%   M:      sparse diagonal matrix, size = [n+1, n+1]
%   Minv:   sparse diagonal matrix, size = [n+1, n+1]

M = sparse(1:n+1, 1:n+1, (1.0)./((1:n+1) - 0.5), n+1, n+1);
Minv = sparse(1:n+1, 1:n+1, (1:n+1) - 0.5, n+1, n+1);

end