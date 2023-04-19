function [val] = Taylor_basis_eval(x, y, p)
% [val] = Taylor_basis_eval(x, y)
%   compute values of Taylor basis up to total degree p at (x,y)
%   should be exact in floating point
% 
% input: 
%   x, y:   evaluation coordinates, scalar
%   p:      total degree, p>=0
% 
% output: 
%   val:    column vector of values, size = [((p+1)*(p+2))/2, 1]

id2D = @(i,j,p) ((2*p+4-i)*(i-1))/2 + j;
N = ((p+1)*(p+2))/2;
val = ones(N, 1);
for i = 1: p+1
    for j = 1: p+2-i
        val(id2D(i,j,p)) = (x^(i-1)) * (y^(j-1));
    end
end

end