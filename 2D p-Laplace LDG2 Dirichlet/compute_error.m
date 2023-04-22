function [total_error] = compute_error(mesh, sol, p, uu, qqx, qqy, basis_type)
% [total_error] = compute_error(mesh, sol, p, uu, qqx, qqy)
%   compute the L2 error of the LDG solution

Np = ((p+1)*(p+2))/2;
[xn, yn, zn] = node_lobatto2D(p, false);
[~, Minv, ~, ~, ~, ~] = matrix_nodal2D(xn, yn, zn, p, basis_type);

% interior quadrature (for the ref element: the right triangle)
[qx, qy, qz, qw, nquad] = quadrature_nodal2D(2*p, false); % algebraic exactness = 2*p
Pnq = interpmat_nodal2D(xn, yn, p, qx, qy, basis_type); % maps nodal values to quadrature nodal values (for quadrature-based interpolation)
Pqn = Minv*Pnq'*diag(qw); % quadrature-based interpolation operator (maps values at quadrature nodes to nodal representations, which is not necessarily the inverse of Pnq)

err_elem = zeros(size(mesh.elems,1),2);  % store error in u and q, respectively
parfor k = 1: size(mesh.elems,1)
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    detJac = det(Jac);

    uval = zeros(nquad,1);
    qxval = zeros(nquad,1);
    qyval = zeros(nquad,1);
    uuval = Pnq * sol.u(:, k);
    qqxval = Pnq * sol.q(1:Np, k);
    qqyval = Pnq * sol.q(Np+1:2*Np, k);
    for i = 1: nquad
        xp = dot(xx, [qx(i), qy(i), qz(i)]);
        yp = dot(yy, [qx(i), qy(i), qz(i)]);
        uval(i) = uu(xp,yp);
        qxval(i) = qqx(xp,yp);
        qyval(i) = qqy(xp,yp);
    end
    err_elem(k,:) = [
        dot(detJac*qw, abs(uuval - uval).^2);
        dot(detJac*qw, abs(qqxval - qxval).^2 + abs(qqyval - qyval).^2)
    ];
end
total_error = [sqrt(sum(err_elem(:,1))), sqrt(sum(err_elem(:,2)))];

end