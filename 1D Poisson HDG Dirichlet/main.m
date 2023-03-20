uu = @(x) exp(sin(x));
qq = @(x) exp(sin(x))*cos(x);
ff = @(x) exp(sin(x))*(sin(x) - cos(x)^2);

a = 0;
b = 3*pi;
gl = uu(a);
gr = uu(b);
type = 'lobatto';   % collocation type
show_spy = false;
tau = 10.0;
post_process = true;

p = 4;

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
    % N(k) = 5*k+10;

    xx = linspace(a, b, N(k)+1);
    hh = xx(2:end) - xx(1:end-1);
    h(k) = max(hh(:));

    [q, u, uhat, qhat] = poisson1D_HDG(ff, xx, gl, gr, p, tau, type, show_spy, 'direct');

    if post_process
        ustar = LDG_postprocess(u, qhat, ff, xx, p, type);
        err(k,:) = zeros(1,4);
        err_loc = zeros(length(xx)-1, 4);
    else
        err(k,:) = zeros(1,3);
        err_loc = zeros(length(xx)-1, 3);
    end

    for j = 1: length(xx) - 1
        inds = (j-1)*(p+1)+1 : j*(p+1);
        indstar = (j-1)*(p+2)+1 : j*(p+2);
        uvals = Pnq*u(inds);
        qvals = Pnq*q(inds);
        ustarvals = Pnqstar*ustar(indstar);
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

        if post_process
            % error of the post-processed solution (should be of order p+2 when p>=1)
            err_loc(j,4) = (hh(j)/2) * dot(qw, abs(uuvals - ustarvals).^2);
        end
    end

    % error in uhat (integrate on each edge only once)
    for j = 1: length(xx)
        err(k,3) = err(k,3) + sqrt(hh(max([j-1,1])) * hh(min([j,length(xx)-1])))*abs(uu(xx(j)) - uhat(j))^2;
    end

    err(k,1) = sum(err_loc(:,1));
    err(k,2) = sum(err_loc(:,2));
    if post_process
        err(k,4) = sum(err_loc(:,4));
    end
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