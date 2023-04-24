% only for p_equ = 2
disp('REMARK: this example is only valid for p_equ == 2, not for general p-Laplace equations!');
uu = @(x,y) exp(sin(pi*x))*cos(pi*(x+y));
qqx = @(x,y) exp(sin(pi*x))*pi*(cos(pi*x)*cos(pi*(x+y))-sin(pi*(x+y)));
qqy = @(x,y) -exp(sin(pi*x))*pi*sin(pi*(x+y));
fluxx = qqx;
fluxy = qqy;
ff = @(x,y) exp(sin(pi*x))*pi*pi*(cos(pi*(x+y))*(2.0-cos(pi*x)^2+sin(pi*x))+2*cos(pi*x)*sin(pi*(x+y)));