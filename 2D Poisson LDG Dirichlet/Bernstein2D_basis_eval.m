function [val, dxval, dyval] = Bernstein2D_basis_eval(x, y, p)
% [val] = Bernstein2D_basis_eval(x, y, p)
%   compute values of triangular Bernstein basis up to total degree p at (x,y)
%   should be accurate in floating point
% 
% input: 
%   x, y:   evaluation coordinates, scalar
%   p:      total degree, p>=0
% 
% output: 
%   val, dxval, dyval:    column vector of values, size = [((p+1)*(p+2))/2, 1]

if p < 0
    error('Bernstein2D_basis_eval: negative p!');
end

if numel(x) ~= 1 || numel(y) ~= 1
    error('Bernstein2D_basis_eval: error x,y shape!');
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
        k = p+3-i-j;
        if nargout == 1
            val(id2D(i,j,p)) = Bernstein2DP(x, y, i-1, j-1, k-1);
        elseif nargout > 1
            [val(id2D(i,j,p)), dxval(id2D(i,j,p)), dyval(id2D(i,j,p))] = Bernstein2DP(x, y, i-1, j-1, k-1);
        end
    end
end

end