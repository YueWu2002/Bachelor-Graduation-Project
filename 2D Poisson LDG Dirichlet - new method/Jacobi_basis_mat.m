function [M, Minv, D] = Jacobi_basis_mat(n, alpha, beta)
% M, Minv, D] = Jacobi_basis_mat(n, alpha, beta)
%   compute the matrices for Jacobi basis on [-1,1]
% 
% input: 
%   n:              max degree, n>=0
%   alpha, beta:    parameters, scalar, > -1
% 
% output: 
%   M:      mass matrix, size = [n+1, n+1]
%   Minv:   inverse of the mass matrix, size = [n+1, n+1]
%   D:      differentiation matrix, size = [n+1, n+1]

if min(alpha, beta) <= -1
    error('Jacobi Polynomial: invalid alpha or beta!');
end

gamma = alpha + beta;
eta = alpha - beta;
zeta = gamma*eta;

M = sparse(1:n+1, 1:n+1, (1.0)./((1:n+1) - 0.5), n+1, n+1, n+1);
Minv = sparse(1:n+1, 1:n+1, (1:n+1) - 0.5, n+1, n+1, n+1);

% p_{2n+1}'(x) = sum_{k=0}^{n-1} (4k+3)*p_{2k+2}(x)
% p_{2n+2}'(x) = sum_{k=0}^{n} (4k+1)*p_{2k+1}(x)
% (index starts from 1)

D = zeros(n+1, n+1);
for j = 1: n+1
    for i = 1: j-1
        if (mod(i+j, 2) ~= 0)
            D(i,j) = 2*i - 1;
        end
    end
end

end