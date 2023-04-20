function [val, dxval, dyval] = Bernstein2DP(x, y, i, j, k)
% Bernstein polynomial on the right triangle with nodes (1,0), (0,1), (0,0)
%   alpha_{i,j,k} * x^(i) * y^(j) * z^(k), z = 1-x-y
% 
% input: 
%   x, y:               arrays of the same shape
%   i, j, k:            integers, scalar
% 
% output: 
%   val, dxval, dyval:  array of the same shape as x,y,z

if (size(x) ~= size(y))
    error('Bernstein2DP: nonconforming array shape!');
end

if (i < 0) || (j < 0) || (k < 0) || nargout < 1
    if nargout == 1
        val = zeros(size(x));
    elseif nargout > 1
        dxval = zeros(size(x));
        dyval = zeros(size(x));
    end
    return;
end

i = double(i);
j = double(j);
k = double(k);
n = i + j + k;

% alpha = gamma(n+1) / (gamma(i+1) * gamma(j+1) * gamma(k+1));
% alpha = double( factorial(int64(n)) / (factorial(int64(i)) * factorial(int64(j)) * factorial(int64(k))) );
alpha = double( nchoosek(int64(n), int64(i)) * nchoosek(int64(j+k), int64(j)) ); % avoid overflow, and increse accuracy (which proves effective in very high-order computation)

z = 1.0 - (x + y);
val = alpha .* ((x.^(i)) .* (y.^(j)) .* (z.^(k)));

if nargout <= 1
    return;
end

% derivatives
dxval = n .* (Bernstein2DP(x, y, i-1, j, k) - Bernstein2DP(x, y, i, j, k-1));
dyval = n .* (Bernstein2DP(x, y, i, j-1, k) - Bernstein2DP(x, y, i, j, k-1));

end