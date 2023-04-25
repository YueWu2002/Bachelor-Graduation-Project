r = @(x,y) sqrt(x^2 + y^2);
uu = @(x,y) r(x,y)^((p_equ-2)/(p_equ-1));
qqx = @(x,y) ((p_equ-2)/(p_equ-1)) * r(x,y)^(-1/(p_equ-1)) * cos(atan2(y,x));
qqy = @(x,y) ((p_equ-2)/(p_equ-1)) * r(x,y)^(-1/(p_equ-1)) * sin(atan2(y,x));
fluxx = @(x,y) sign(p_equ-2) * abs((p_equ-2)/(p_equ-1))^(p_equ-1) * (r(x,y))^(-1) * cos(atan2(y,x));
fluxy = @(x,y) sign(p_equ-2) * abs((p_equ-2)/(p_equ-1))^(p_equ-1) * (r(x,y))^(-1) * sin(atan2(y,x));
ff = @(x,y) 0.0;