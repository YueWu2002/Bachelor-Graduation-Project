% reference: Finite element approximation of the p-Laplacian
r = @(x,y) sqrt(x^2 + y^2);
uu = @(x,y) ((p_equ-1)/((sigma+2)^(1/(p_equ-1)))) * (1-r(x,y)^((sigma+p_equ)/(p_equ-1)))/(sigma+p_equ);
qqx = @(x,y) -r(x,y)^((sigma+1)/(p_equ-1)) / (sigma+2)^(1/(p_equ-1)) * cos(atan2(y,x));
qqy = @(x,y) -r(x,y)^((sigma+1)/(p_equ-1)) / (sigma+2)^(1/(p_equ-1)) * sin(atan2(y,x));
fluxx = @(x,y) -r(x,y)^sigma / (sigma+2) * x;
fluxy = @(x,y) -r(x,y)^sigma / (sigma+2) * y;
ff = @(x,y) r(x,y)^sigma;
