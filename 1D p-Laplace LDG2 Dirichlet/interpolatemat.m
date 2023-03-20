function [P] = interpolatemat(x0, x1)
% [P] = interpolatemat(x0, x1)
%   compute the Lagrange polynomial interpolation matrix
% 
% input: 
%   x0: coordinates of the original data set
%   x1: coordinates of the target data set
% 
% output: 
%   P:  interpolation matrix, size = [numel(x1), numel(x0)]

m = numel(x1);
n = numel(x0);
P = ones(m,n);
x0 = reshape(x0, [n,1]);
x1 = reshape(x1, [m,1]);

for j = 1: n
    for k = 1: n
        if (k~=j) 
            P(:,j) = P(:,j).*((x1(:) - x0(k))/(x0(j) - x0(k)));
        end
    end
end

end