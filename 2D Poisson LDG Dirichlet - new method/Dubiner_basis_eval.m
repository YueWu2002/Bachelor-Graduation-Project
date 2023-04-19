function [val] = Dubiner_basis_eval(x, y, p)
% [val] = Dubiner_basis_eval(x, y, p)
%   compute values of Dubiner basis up to total degree p at (x,y)
%   should be accurate in floating point
% 
% input: 
%   x, y:   evaluation coordinates, scalar
%   p:      total degree, p>=0
% 
% output: 
%   val:    column vector of values, size = [((p+1)*(p+2))/2, 1]

id2D = @(i,j,p) ((2*p+4-i)*(i-1))/2 + j; % sequence: i from 1 to p+1, j from 1 to p+2-i
N = ((p+1)*(p+2))/2;
val = ones(N, 1);



end