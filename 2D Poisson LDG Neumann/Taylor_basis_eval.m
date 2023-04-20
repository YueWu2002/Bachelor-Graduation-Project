function [val, dxval, dyval] = Taylor_basis_eval(x, y, p)
% [val] = Taylor_basis_eval(x, y, p)
%   compute values of Taylor basis up to total degree p at (x,y)
%   should be exact in floating point
% 
% input: 
%   x, y:   evaluation coordinates, scalar
%   p:      total degree, p>=0
% 
% output: 
%   val, dxval, dyval:    column vector of values, size = [((p+1)*(p+2))/2, 1]

if p < 0
    error('Taylor_basis_eval: negative p!');
end

if numel(x) ~= 1 || numel(y) ~= 1
    error('Taylor_basis_eval: error x,y shape!');
end

id2D = @(i,j,p) ((2*p+4-i)*(i-1))/2 + j; % sequence: i from 1 to p+1, j from 1 to p+2-i
N = ((p+1)*(p+2))/2;

if nargout == 1
    val = zeros(N, 1);
elseif nargout > 1
    dxval = zeros(N, 1);
    dyval = zeros(N, 1);
end

for i = 1: p+1
    for j = 1: p+2-i
        if nargout == 1
            val(id2D(i,j,p)) = x^(i-1) * y^(j-1);
        elseif nargout > 1
            if i == 1
                dxval(id2D(i,j,p)) = 0.0;
            else
                dxval(id2D(i,j,p)) = (i-1) * x^(i-2) * y^(j-1);
            end
            if j == 1
                dyval(id2D(i,j,p)) = 0.0;
            else
                dyval(id2D(i,j,p)) = (j-1) * x^(i-1) * y^(j-2);
            end
        end
    end
end

end