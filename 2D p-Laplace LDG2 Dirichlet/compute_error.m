function [total_error] = compute_error(mesh, sol, p, p_equ, uu, qqx, qqy, fluxx, fluxy, basis_type)
% [total_error] = compute_error(mesh, sol, p, p_equ, uu, qqx, qqy, fluxx, fluxy, basis_type)
%   compute the error of the LDG solution
%   Lp for u, Lp for q, Lq for sigma

q_equ = p_equ / (p_equ-1.0);

Np = ((p+1)*(p+2))/2;

[xn, yn, zn] = node_lobatto2D(p, false);
[~, Minv, ~, ~, ~, ~] = matrix_nodal2D(xn, yn, zn, p, basis_type);

% interior quadrature (for the ref element: the right triangle)
[qx, qy, qz, qw, nquad] = quadrature_nodal2D(2*p, false); % algebraic exactness = 2*p
Pnq = interpmat_nodal2D(xn, yn, p, qx, qy, basis_type); % maps nodal values to quadrature nodal values (for quadrature-based interpolation)
Pqn = Minv*Pnq'*diag(qw); % quadrature-based interpolation operator (maps values at quadrature nodes to nodal representations, which is not necessarily the inverse of Pnq)

err_elem = zeros(size(mesh.elems,1), 3);  % store error in (u, q, flux)
parfor k = 1: size(mesh.elems,1)
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    detJac = det(Jac);

    % numerical solution at the quadrature nodes
    uuval = Pnq * sol.u(:,k);
    qqval = Pnq * sol.q(:,:,k);
    ssval = Pnq * sol.flux(:,:,k);

    % exact solution
    uval = nan(nquad, 1);
    qval = nan(nquad, 2);
    sval = nan(nquad, 2);
    for i = 1: nquad
        xp = dot(xx, [qx(i), qy(i), qz(i)]);
        yp = dot(yy, [qx(i), qy(i), qz(i)]);
        uval(i) = uu(xp,yp);
        qval(i,:) = [qqx(xp,yp), qqy(xp,yp)];
        sval(i,:) = [fluxx(xp,yp), fluxy(xp,yp)];
    end
    err_elem(k,:) = [
        dot(detJac*qw, abs(uuval(:) - uval(:)).^p_equ); 
        dot(detJac*qw, (abs(qqval(:,1) - qval(:,1)).^2 + abs(qqval(:,2) - qval(:,2)).^2).^(p_equ/2.0)); 
        dot(detJac*qw, (abs(ssval(:,1) - sval(:,1)).^2 + abs(ssval(:,2) - sval(:,2)).^2).^(q_equ/2.0))
    ];
end
total_error = [sum(err_elem(:,1))^(1.0/p_equ), sum(err_elem(:,2))^(1.0/p_equ), sum(err_elem(:,3))^(1.0/q_equ)];

end