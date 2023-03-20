function [q, u] = poisson1D_LDG(f, xx, gl, gr, c0, p, alpha, C12, type, show_spy, solve_method)
% 1D poisson equation with purely Neumann BC
% p: polynomial order
% alpha: C_{11} penalty order, 1 or 0

[qx, qw] = quadrature_nodal1D(p+1, 'legendre');     % quadrature for compute source integration
[nx, ~] = quadrature_nodal1D(p+1, type);            % collocation
[M, Minv, D, ~, ~, Tl, Tr] = matrix_nodal1D(nx);    % collocation
Pnq = interpolatemat(nx, qx);   % maps nodal values at collocation nodes to those at quadrature nodes

N = length(xx) - 1; % number of elems
hh = xx(2:end) - xx(1:end-1);

avg = [0.5*Tr, 0.5*Tl];
jump = [Tr, -Tl];
ETA = 1.0*ones(N+1, 1);


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
    B_block{k} = sparse(-((hh(k)/2) *M) * ((2/hh(k)) * D));
end
B = blkdiag(B_block{:});
% add jump and avg terms (B[i,j] block is B[u_j, tau_i]) (loop over interior face)
for i = 1: N-1    % face between I_{i} and I_{i+1}
    ind = (i-1)*(p+1)+1 : (i+1)*(p+1);
    B(ind, ind) = B(ind, ind) + (avg - C12*jump)'*jump;
end
% no Dirichlet boundary terms


% assemble C
C = sparse(N*(p+1), N*(p+1));
% add jump terms
for i = 1: N-1
    % between I_{i} and I_{i+1}
    ind = (i-1)*(p+1)+1 : (i+1)*(p+1);
    C11 = ETA(i+1) / sqrt(hh(i) * hh(i+1))^alpha;
    C(ind, ind) = C(ind, ind) + jump'*C11*jump;
end
% no Dirichlet boundary terms


% assemble F (only contains Dirichlet data)
F = sparse(N*(p+1), 1);
% no Dirichlet boundary terms


% assemble G
G = zeros(N*(p+1), 1);
for k = 1: N
    base = (k-1)*(p+1);
    xc = 0.5*(xx(k) + xx(k+1));

    % quadrature-based interpolation
    fvals = zeros(length(qx),1);
    for j = 1: length(qx)
        xp = xc + (hh(k)/2)*qx(j);
        fvals(j) = f(xp);
    end
    G(base+1 : base+(p+1)) = Pnq'*((hh(k)/2)*diag(qw))*fvals;
end
% no Dirichlet boundary terms
% add Neumann boundary terms
ind = 1: p+1;
G(ind) = G(ind) + gl*Tl';
ind = (N-1)*(p+1)+1 : N*(p+1);
G(ind) = G(ind) + gr*Tr';


% assemble  Lm (the Lagrange multiplier term)
Lm = A * ones(N*(p+1), 1);


% solve
switch solve_method
case 'direct'
    Z = [
        A, B, sparse(N*(p+1), 1);
        B', -C, Lm
        sparse(1, N*(p+1)), Lm', 0
    ];
    x = Z \ [F; -G; c0];
    q = x(1:N*(p+1));
    u = x(N*(p+1)+1 : 2*N*(p+1));
    if show_spy
        spy(Z);
        print(gcf, '-depsc2', '-vector', 'spyZZ', '-r300');
    end
case 'reduce'
    S = [
        C + B'*Ainv*B, Lm;
        Lm', 0
    ];
    s1 = Ainv*F;
    s2 = S \ [B'*s1 + G; c0];
    u = s2(1:N*(p+1));
    q = s1 - Ainv*(B*u);
    if show_spy
        spy(S);
        print(gcf, '-depsc2', '-vector', 'spySS', '-r300');
    end
end

end