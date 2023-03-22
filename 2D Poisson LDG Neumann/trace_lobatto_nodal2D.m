function [T] = trace_lobatto_nodal2D(p)
% [T] = trace_lobatto_nodal2D(p)
%   maps the interior nodal DOFs at triangular Lobatto points to boundary nodal DOFs at edge Lobatto points
%   domain: the right triangle
%   edge sequence: 
%       EDGE1: (0,1)->(0,0);
%       EDGE2: (0,0)->(1,0);
%       EDGE3: (1,0)->(0,1);
% 
% input: 
%   p:  the polynomial total degree
% 
% output: 
%   T1, T2, T3: mapping matrices, size = [p+1, ((p+1)*(p+2))/2] (we combine them into a dense tensor)
%
% REMARK: please make the matrix sparse before putting it into computing!
% WARNING: please plot the nodal distribution WITH ORDER NUMBERS before writing such code!

% checked.

id2D = @(i,j,p) ((2*p+4-i)*(i-1))/2 + j;
N = ((p+1)*(p+2))/2;
T = zeros(p+1, N, 3);

for k = 1: p+1
    T(p+2-k, k, 1) = 1;
    T(k, id2D(k,1,p), 2) = 1;
    T(p+2-k, id2D(k,p+2-k,p), 3) = 1;
end

end