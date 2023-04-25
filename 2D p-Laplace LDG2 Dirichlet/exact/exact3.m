
uu = @(x,y) u(x,y,a);
qqx = @(x,y) q(x,y,a) * cos(atan2(y,x));
qqy = @(x,y) q(x,y,a) * sin(atan2(y,x));
fluxx = @(x,y) q(x,y,a)^(p_equ-1) * cos(atan2(y,x));
fluxy = @(x,y) q(x,y,a)^(p_equ-1) * sin(atan2(y,x));
ff = @(x,y) f(x,y,a,p_equ);

function [val] = r(x,y)
val = sqrt(x^2 + y^2);
end

function [val] = u(x,y,a)
if r(x,y) > a
    val = (r(x,y)-a)^4;
else
    val = 0.0;
end
end

function [val] = q(x,y,a)
if r(x,y) > a
    val = 4*(r(x,y)-a)^3;
else
    val = 0.0;
end
end

function [val] = f(x,y,a, p_equ)
if r(x,y) > a
    val = 4^(p_equ-1) * (r(x,y)-a)^(3*p_equ-4) * (2-3*p_equ+a/r(x,y));
else
    val = 0.0;
end
end