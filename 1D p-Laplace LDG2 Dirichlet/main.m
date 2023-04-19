p_equ = 1.3;
q_equ = p_equ/(p_equ-1);

uu = @(x) (x + 1.0).^((p_equ-2)/(p_equ-1)); % u(x)
dduu = @(x) ((p_equ-2)/(p_equ-1)) * (x + 1.0).^(-1.0/(p_equ - 1));    % u'(x)
qq = @(x) (abs((p_equ-2)/(p_equ-1))^(p_equ-1))*sign(p_equ-2)./(x+1);    % flux[u](x)
ff = @(x) (abs((p_equ-2)/(p_equ-1))^(p_equ-1))*sign(p_equ-2)./((x+1).^2);  % f(x)

a = 0;
b = 3*pi;
gl = uu(a);
gr = uu(b);
type = 'lobatto';   % collocation type (should have no impact on the solution)
show_spy = false;

p = 4;

[qx, qw] = quadrature_nodal1D(p+1, 'legendre');     % quadrature for detect error
[nx, ~] = quadrature_nodal1D(p+1, type);            % collocation
[M, Minv, D, ~, ~, Tl, Tr] = matrix_nodal1D(nx);    % collocation
Pnq = interpolatemat(nx, qx);                       % maps nodal values at collocation nodes to those at quadrature nodes
Pqn = Minv*Pnq'*diag(qw);                           % quadrature-based interpolation operator (maps values at quadrature nodes to nodal representations)

err = [];
h = [];
N = [];
for k = 1: 6
     N(k) = 2^(k+3);
    % N(k) = 8+8*k;
    % N(k) =  16+4*k;

    xx = linspace(a, b, N(k)+1);
    hh = xx(2:end) - xx(1:end-1);
    h(k) = max(hh(:));

    [u, du, q] = pLaplace1D_LDG(ff, xx, gl, gr, p, p_equ, 1, 0.0, type, true);
    % results: nodal representation of the displacement u and its gradient du
    
    err(k,:) = [0, 0, 0];
    err_loc = zeros(length(xx)-1, 3);
    for j = 1: length(xx) - 1
        inds = (j-1)*(p+1)+1 : j*(p+1);
        uvals = Pnq*u(inds);                % evaluate polynomial at quadrature nodes
        duvals = Pnq*du(inds);              % evaluate polynomial at quadrature nodes
        % qvals = Pnq*(Pqn*Amap(duvals, p_equ)); % quadrature-based interpolation + evaluate at quadrature points
        qvals = Pnq*q(inds);

        uuvals = zeros(length(qw),1);       % exact solution
        dduuvals = zeros(length(qw),1);     % exact solution
        qqvals = zeros(length(qw),1);       % exact solution

        xc = 0.5*(xx(j) + xx(j+1));
        for i = 1: length(qx)
            xp = xc + 0.5*hh(j)*qx(i);
            uuvals(i) = uu(xp);
            dduuvals(i) = dduu(xp);
            qqvals(i) = qq(xp);
        end
        err_loc(j,1) = dot((hh(j)/2)*qw, abs(uvals - uuvals).^p_equ);       % Lp
        err_loc(j,2) = dot((hh(j)/2)*qw, abs(duvals - dduuvals).^p_equ);    % Lp
        err_loc(j,3) = dot((hh(j)/2)*qw, abs(qvals - qqvals).^q_equ);       % Lq
    end

    err(k,1) = sum(err_loc(:,1))^(1.0/p_equ);
    err(k,2) = sum(err_loc(:,2))^(1.0/p_equ);
    err(k,3) = sum(err_loc(:,3))^(1.0/q_equ);
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