function [y] = Amap(x, p_equ)
% y = x^(p-1)*sgn(x), if p~=1
y = sign(x).*(abs(x).^(p_equ-1));
end