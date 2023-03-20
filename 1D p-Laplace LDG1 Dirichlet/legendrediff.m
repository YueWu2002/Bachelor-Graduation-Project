function [D] = legendrediff(n)
% [D] = legendrediff(n)
%   compute the differentiation matrix with respect to Legendre polynomials
%
% input: 
%   n:  max order of Legendre polynomials
% 
% output: 
%   D:      strictly upper-triangular matrix, size = [n+1, n+1], D(i,j): coefficient of p(i) in p'(j)

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