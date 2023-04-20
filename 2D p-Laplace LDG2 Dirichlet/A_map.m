function [y] = A_map(x, p_equ)
% [y] = A_map(x, p_equ)
%   y = |x|^(p-1) * e(x), where x is a vector, e(x) is the direction of x
% 
% input: 
%   x:      column or row vector
%   p_equ:  scalar, p_equ >= 1
% 
% output: 
%   y:      array, the same shape as x

if p_equ <= 1.0
    error('A_map: error p!');
end

if ~isvector(x)
    error('A_map: error shape of x!');
end

rho = norm(x, 2);
rho_ans = rho^(p_equ-1);
if rho_ans > 0.0
    x_dir = x / rho;
    y = rho_ans * x_dir;
else
    y = zeros(size(x));
end

end