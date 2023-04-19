function [y] = jacobip(n, alpha, beta, x)
% [y] = jacobip(n, alpha, beta, x)
% 
% input: 
%   n:              order, scalar, n>=0
%   alpha, beta:    parameters, scalar, > -1
%   x:              array
% 
% output: 
%   y:              values, the same shape as x
% 
% References: 
%   [1] https://arxiv.org/abs/1509.08213
%   [2] Shen, J., Tang, T., Wang, LL. (2011). Orthogonal Polynomials and Related Approximation Results. In: Spectral Methods. Springer Series in Computational Mathematics, vol 41. Springer, Berlin, Heidelberg. https://doi.org/10.1007/978-3-540-71041-7_3
%   [3] https://people.sc.fsu.edu/~jburkardt/f77_src/jacobi_polynomial/jacobi_polynomial.html

if min(alpha, beta) <= -1
    error('Jacobi Polynomial: invalid alpha or beta!');
end

gamma = alpha + beta;
eta = alpha - beta;
zeta = gamma*eta;

if n < 0
    y = zeros(size(x));
    return;
elseif n == 0
    y = ones(size(x));
    return;
elseif n == 1
    y = (1.0 + 0.5*gamma) .* x + 0.5*eta;
    return;
else
    t1 = zeros(size(x)); % k=-1
    t2 = ones(size(x)); % k=0
    y = (1.0 + 0.5*gamma) * x + 0.5*eta; % k=1
    for k = 2: n
        t1 = t2;
        t2 = y;

        % k >= 1
        if k == 1
            A = (2*k + gamma) / (2*k);
            B = eta / (2*k + gamma);
        else
            A = (2*k + gamma) * (2*k - 1  + gamma) / (2*k*(k + gamma));
            B = zeta / ((2*k + gamma - 2)*(2*k + gamma));
        end
        % to be continued
        C = 2*(k + alpha - 1)*(k + beta - 1) / ((2*k + gamma - 2)*(2*k + gamma - 1));
        y = A .* ((x + B) .* t2 - C .* t1);
        % P_{n} = A * ((x + B)*P_{n-1} - C * P_{n-2})
    end
end

end