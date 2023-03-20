uu = @(x) exp(sin(x));
qq = @(x) exp(sin(x))*cos(x);
ff = @(x) exp(sin(x))*(sin(x) - cos(x)^2);

a = 0;
b = 3*pi;
c0 = 14.16368455672395547087675939170498;
gl = -qq(a);    % directional derivative
gr = qq(b);
type = 'lobatto';
show_spy = false;

p = 2;

[qx, qw] = quadrature_nodal1D(p+1, 'legendre');     % quadrature for detect error
[nx, ~] = quadrature_nodal1D(p+1, type);            % collocation
[M, Minv, D, ~, ~, Tl, Tr] = matrix_nodal1D(nx);    % collocation
Pnq = interpolatemat(nx, qx);                       % maps nodal values at collocation nodes to those at quadrature nodes

[nxstar, ~] = quadrature_nodal1D(p+2, type);        % collocation (for post-processed solution)
Pnqstar = interpolatemat(nxstar, qx);               % maps nodal values at collocation nodes to those at quadrature nodes

err = [];
h = [];
N = [];
for k = 1: 6
    N(k) = 2^(k+3);

    xx = linspace(a, b, N(k)+1);
    hh = xx(2:end) - xx(1:end-1);
    h(k) = max(hh(:));

    [q, u] = poisson1D_LDG(ff, xx, gl, gr, c0, p, 1, 0.0, type, show_spy, 'direct');
    
    err(k,:) = zeros(1,2);
    err_loc = zeros(length(xx)-1, 2);
    for j = 1: length(xx) - 1
        inds = (j-1)*(p+1)+1 : j*(p+1);
        indstar = (j-1)*(p+2)+1 : j*(p+2);
        uvals = Pnq*u(inds);
        qvals = Pnq*q(inds);
        uuvals = zeros(length(qw),1);   % exact solution
        qqvals = zeros(length(qw),1);   % exact solution

        xc = 0.5*(xx(j) + xx(j+1));
        for i = 1: length(qx)
            xp = xc + 0.5*hh(j)*qx(i);
            uuvals(i) = uu(xp);
            qqvals(i) = qq(xp);
        end
        err_loc(j,1) = dot((hh(j)/2)*qw, abs(uvals - uuvals).^2);
        err_loc(j,2) = dot((hh(j)/2)*qw, abs(qvals - qqvals).^2);
    end

    err(k,1) = sum(err_loc(:,1));
    err(k,2) = sum(err_loc(:,2));

    err(k,:) = sqrt(err(k,:));
end

ord = nan(size(err));
for k = 2: size(ord, 1)
    ord(k, :) = log(err(k,:)./err(k-1,:))/log(h(k)/h(k-1));
end

%{
format short e
disp(err)

format short
disp(ord)
%}

for i = 1: size(err, 1)
    fprintf('%d ', N(i));
    for j = 1: size(err, 2)
        fprintf("& %.2e & %.2f ", [err(i,j), ord(i,j)]);
    end
    fprintf('\\\\\n');
end