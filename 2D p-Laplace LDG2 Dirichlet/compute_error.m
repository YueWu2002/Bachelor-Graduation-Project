function [total_error] = compute_error(mesh, sol, p, uu, qqx, qqy)
% [total_error] = compute_error(mesh, sol, p, uu, qqx, qqy)
%   compute the L2 error of the LDG solution

Np = ((p+1)*(p+2))/2;
[xn, yn, zn] = node_lobatto2D(p, false);
[M, ~, ~, ~, Pmn, ~] = matrix_nodal2D(xn, yn, zn, p);

% interior quadrature (for the ref element: the right triangle)
[qx, qy, qz, qw, nquad] = quadrature_nodal2D(2*p, false); % algebraic exactness = 2*p
Pmq = zeros(nquad, Np); % maps Taylor basis coefficients to quadrature nodal values
for k = 1: nquad
    Pmq(k,:) = Taylor_basis_eval(qx(k), qy(k), p);
end
Pnq = Pmq / Pmn;  % maps nodal values to quadrature nodal values (for quadrature-based interpolation)
% Pnq = Pmq*Pnm;

err_elem = zeros(size(mesh.elems,1),2);  % store error in u and q, respectively
parfor k = 1: size(mesh.elems,1)
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    M_elem = det(Jac)*M;

    %{
    % interpolation-based quadrature
    ueval = zeros(Np, 1);
    qeval = zeros(Np, 2);
    for i = 1: Np
        xp = dot(xx, [xn(i), yn(i), zn(i)]);
        yp = dot(yy, [xn(i), yn(i), zn(i)]);
        ueval(i) = uu(xp,yp) - sol.u((k-1)*Np + i);
        qeval(i,:) = [qqx(xp,yp) - sol.q((k-1)*(2*Np)+i), qqy(xp,yp) - sol.q((k-1)*(2*Np)+i+Np)];
    end
    err_elem(k,:) = [ueval'*M_elem*ueval, qeval(:,1)'*M_elem*qeval(:,1) + qeval(:,2)'*M_elem*qeval(:,2)];
    %}

    % quadrature
    uval = zeros(nquad,1);
    qxval = zeros(nquad,1);
    qyval = zeros(nquad,1);
    uuval = Pnq * sol.u((k-1)*Np+1 : k*Np);
    qqxval = Pnq * sol.q((k-1)*(2*Np)+1 : (k-1)*(2*Np)+Np);
    qqyval = Pnq * sol.q((k-1)*(2*Np)+Np+1 : k*(2*Np));
    for i = 1: nquad
        xp = dot(xx, [qx(i), qy(i), qz(i)]);
        yp = dot(yy, [qx(i), qy(i), qz(i)]);
        uval(i) = uu(xp,yp);
        qxval(i) = qqx(xp,yp);
        qyval(i) = qqy(xp,yp);
    end
    err_elem(k,:) = [
        dot(det(Jac)*qw, abs(uuval - uval).^2);
        dot(det(Jac)*qw, abs(qqxval - qxval).^2 + abs(qqyval - qyval).^2)
    ];
end
total_error = [sqrt(sum(err_elem(:,1))), sqrt(sum(err_elem(:,2)))];

end