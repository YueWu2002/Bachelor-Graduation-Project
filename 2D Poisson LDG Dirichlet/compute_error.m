function [total_error] = compute_error(mesh, sol, p, uu, qqx, qqy)
% [total_error] = compute_error(mesh, sol, p, uu, qqx, qqy)
%   compute the L2 error of the LDG solution

[xn, yn, zn] = node_lobatto2D(p, false);
[M, ~, ~, ~, ~, ~] = matrix_nodal2D(xn, yn, zn, p);
Np = ((p+1)*(p+2))/2;
id2D = @(i,j,s) ((2*s+4-i)*(i-1))/2 + j;
err_elem = zeros(size(mesh.elems,1),2);  % store error in u and q, respectively
parfor k = 1: size(mesh.elems,1)
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    M_elem = det(Jac)*M;
    ueval = zeros(Np, 1);
    qeval = zeros(Np, 2);
    for i = 1: Np
        xp = dot(xx, [xn(i), yn(i), zn(i)]);
        yp = dot(yy, [xn(i), yn(i), zn(i)]);
        ueval(i) = uu(xp,yp) - sol.u((k-1)*Np + i);
        qeval(i,:) = [qqx(xp,yp) - sol.q((k-1)*(2*Np)+i), qqy(xp,yp) - sol.q((k-1)*(2*Np)+i+Np)];
    end
    err_elem(k,:) = [ueval'*M_elem*ueval, qeval(:,1)'*M_elem*qeval(:,1) + qeval(:,2)'*M_elem*qeval(:,2)];
end
total_error = [sqrt(sum(err_elem(:,1))), sqrt(sum(err_elem(:,2)))];

end