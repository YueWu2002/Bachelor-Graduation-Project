function [q, u, uhat, qhat] = poisson1D_HDG(f, xx, gl, gr, p, tau, type, show_spy, solve_method)
% 1D poisson equation with purely Dirichlet BC
% p: polynomial order

[qx, qw] = quadrature_nodal1D(p+1, 'legendre');     % quadrature for compute source integration
[nx, ~] = quadrature_nodal1D(p+1, type);            % collocation
[M, Minv, D, ~, ~, Tl, Tr] = matrix_nodal1D(nx);    % collocation
Pnq = interpolatemat(nx, qx);   % maps nodal values at collocation nodes to those at quadrature nodes

N = length(xx) - 1; % number of elems
hh = xx(2:end) - xx(1:end-1);

avg = [0.5*Tr, 0.5*Tl];
jump = [Tr, -Tl];


% assemble Ainv
Ainv_block = cell(N,1);
for k = 1: N
    Ainv_block{k} = sparse((2/hh(k)) * Minv);
end
Ainv = blkdiag(Ainv_block{1:N});


% assemble A
A_block = cell(N,1);
for k = 1: N
    A_block{k} = sparse((hh(k)/2) * M);
end
A = blkdiag(A_block{1:N});


% assemble B
B_block = cell(N,1);
for k = 1: N
    B_block{k} = sparse(((2/hh(k)) * D)' * ((hh(k)/2) * M));
end
B = blkdiag(B_block{1:N});


% assemble C
C_block = cell(N,1);
for k = 1: N
    C_block{k} = sparse((Tl'*tau*Tl) + (Tr'*tau*Tr));
end
C = blkdiag(C_block{1:N});


% assemble X
X = sparse(N+1, N*(p+1));
for k = 1: N
    ind = (k-1)*(p+1)+1 : k*(p+1);
    X(k, ind) = -Tl;
    X(k+1, ind) = Tr;
end


%assemble Y
Y = sparse(N*(p+1), N+1);
for k = 1: N
    ind = (k-1)*(p+1)+1 : k*(p+1);
    Y(ind, k) = tau*Tl;
    Y(ind, k+1) = tau*Tr;
end


% assemble Z
Z = blkdiag(tau, 2*tau*speye(N-1), tau);
% here, pay attention to the '2', because every interior edge was computed twice!!!


% assemble F (source vector as an dual element, use quadrature-based interpolation)
F = zeros(N*(p+1), 1);
for k = 1: N
    base = (k-1)*(p+1);
    xc = 0.5*(xx(k) + xx(k+1));

    % quadrature-based interpolation
    fvals = zeros(length(qx),1);
    for j = 1: length(qx)
        xp = xc + (hh(k)/2)*qx(j);
        fvals(j) = f(xp);
    end
    F(base+1 : base+(p+1)) = Pnq'*((hh(k)/2)*diag(qw))*fvals;
end


% assemble G (only contains Neumann data)
G = sparse(N+1, 1);


% assemble the big matrix equation
ZZ = [
    A, B, -X';
    -B', C, -Y;
    X, -Y', Z
];
RHS = [sparse(N*(p+1), 1); F; G];


switch solve_method
case 'direct'
    % enforce the Dirichlet boundary condition by splitting the matrix

    % id for Dirichlet DOFs
    id_bc = [2*N*(p+1) + 1, 2*N*(p+1) + N+1];

    % id for other DOFs
    id_in = 1 : 2*N*(p+1) + N+1;
    id_in(id_bc) = [];

    % split the equation
    x = ZZ(id_in, id_in) \ (RHS(id_in) - ZZ(id_in, id_bc)*[gl; gr]);
    q = x(1:N*(p+1));
    u = x(N*(p+1)+1 : 2*N*(p+1));
    uhat = [gl; x(2*N*(p+1)+1 : end); gr];

    if show_spy
        spy(ZZ);
        print(gcf, '-depsc2', '-vector', 'spyZZ', '-r300');
    end

case 'reduce'
    % assemble TT = C + B'*A^{-1}*B (block diagonal, symmetric positive definite)
    TTinv_block = cell(N,1);
    for k = 1: N
        TTinv_block{k} = inv(C_block{k} + B_block{k}'*Ainv_block{k}*B_block{k});
    end
    TTinv = blkdiag(TTinv_block{1:N});

    % assemble SS = A^{-1}*B (block structure)
    SS_block = cell(N,1);
    for k = 1: N
        SS_block{k} = Ainv_block{k}*B_block{k};
    end
    SS = blkdiag(SS_block{1:N});

    % assemble matrix for uhat
    L = Y' + X * SS;
    KK = Z + X * Ainv * X' - L * TTinv * L';
    FF = G + L * (TTinv * F);


    % solve uhat with Dirichlet BC
    % id for Dirichlet DOFs and other DOFs
    id_bc = [1, N+1];
    id_in = 1 : N+1;
    id_in(id_bc) = [];

    uhat = [gl; KK(id_in, id_in) \ (FF(id_in) - KK(id_in, id_bc)*[gl; gr]); gr];
    u = TTinv * (L' * uhat + F);
    q = Ainv * (X'*uhat - B*u);

    if show_spy
        spy(KK);
        print(gcf, '-depsc2', '-vector', 'spySS', '-r300');
    end
end

% compute the numerical trace for q
qhat = zeros(N+1, 1);
ind = 1 : p+1;
qhat(1) = Tl*q(ind) - tau*(-1)*(Tl*u(ind) - uhat(1));
ind = (N-1)*(p+1)+1 : N*(p+1);
qhat(N+1) = Tr*q(ind) - tau*(Tr*u(ind) - uhat(N+1));
for k = 1: N-1
    inds = (k-1)*(p+1)+1 : (k+1)*(p+1);
    qhat(k+1) = avg*q(inds) - 0.5*tau*jump*u(inds);
end


end