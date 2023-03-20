function [ustar] = LDG_postprocess(u, qhat, f, xx, p, type)
% post-process the LDG or LDG-H solution
% from [u] in P^{p} to [ustar] in P^{p+1}

[qx, qw] = quadrature_nodal1D(p+1, 'legendre');     % quadrature for compute source integration
[nx, ~] = quadrature_nodal1D(p+1, type);            % collocation
[M, ~, ~, ~, ~, ~, ~] = matrix_nodal1D(nx);       % collocation

[nxstar, ~] = quadrature_nodal1D(p+2, type);        % collocation (for post-process)
[Mstar, ~, Dstar, ~, ~, Tlstar, Trstar] = matrix_nodal1D(nxstar);    % collocation (for post-process)
Pnqstar = interpolatemat(nxstar, qx);               % maps nodal values at collocation nodes to those at quadrature nodes

N = length(xx) - 1;
hh = xx(2:end) - xx(1:end-1);
ustar = zeros(N*(p+2), 1);

for k = 1: N
    ind = (k-1)*(p+1)+1 : k*(p+1);
    indstar = (k-1)*(p+2)+1 : k*(p+2);

    % assemble F
    % quadrature-based interpolation
    fvals = zeros(length(qx),1);
    xc = 0.5*(xx(k) + xx(k+1));
    for j = 1: length(qx)
        xp = xc + (hh(k)/2)*qx(j);
        fvals(j) = f(xp);
    end
    F = Pnqstar'*((hh(k)/2)*diag(qw))*fvals;

    % assemble G (Neumann boundary)
    G = Tlstar' * (-1)*qhat(k) + Trstar' * qhat(k+1);

    % assemble Lm (Lagrange multiplier)
    Lm = ((hh(k)/2)*Mstar) * ones(p+2, 1);

    % assemble S (stiffness matrix)
    S = ((2/hh(k)) * Dstar)' * ((hh(k)/2) * Mstar) * ((2/hh(k)) * Dstar);

    % assemble the integration of u
    c0 = ones(1, p+1)*((hh(k)/2)*M)*u(ind);

    temp = [S, Lm; Lm', 0] \ [F + G; c0];
    ustar(indstar) = temp(1:p+2);

end


end