function [q, u] = poisson1D_LDG(f, xx, gl, gr, p, alpha, C12, type, show_spy, solve_method)
% 1D poisson equation with purely Dirichlet BC
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
% add Dirichlet boundary terms
ind = 1: p+1;
B(ind, ind) = B(ind, ind) - Tl'*Tl;
ind = (N-1)*(p+1)+1 : N*(p+1);
B(ind, ind) = B(ind, ind) + Tr'*Tr;


% assemble C
C = sparse(N*(p+1), N*(p+1));
% add jump terms
for i = 1: N-1
    % between I_{i} and I_{i+1}
    ind = (i-1)*(p+1)+1 : (i+1)*(p+1);
    C11 = ETA(i+1) / sqrt(hh(i) * hh(i+1))^alpha;
    C(ind, ind) = C(ind, ind) + jump'*C11*jump;
end
% add Dirichlet boundary penalty terms
C11 = ETA(1) / hh(1)^alpha;
ind = 1: p+1;
C(ind, ind) = C(ind, ind) + Tl'*C11*Tl;
C11 = ETA(N+1) / hh(N)^alpha;
ind = (N-1)*(p+1)+1 : N*(p+1);
C(ind, ind) = C(ind, ind) + Tr'*C11*Tr;


% assemble F (only contains Dirichlet data)
F = sparse(N*(p+1), 1);
ind = 1: p+1;
F(ind) = F(ind) - gl*Tl';
ind = (N-1)*(p+1)+1 : N*(p+1);
F(ind) = F(ind) + gr*Tr';


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
% add Dirichlet boundary penalty terms
C11 = ETA(1) / hh(1)^alpha;
ind = 1: p+1;
G(ind) = G(ind) + C11*gl*Tl';
C11 = ETA(N+1) / hh(N)^alpha;
ind = (N-1)*(p+1)+1 : N*(p+1);
G(ind) = G(ind) + C11*gr*Tr';
% no Neumann boundary terms


% solve
switch solve_method
case 'direct'
    Z = [A, B; -B', C];
    x = Z \ [F; G];
    q = x(1:N*(p+1));
    u = x(N*(p+1)+1 : end);
    if show_spy
        spy(Z);
        print(gcf, '-depsc2', '-vector', 'spyZZ', '-r300');
    end
case 'reduce'
    % u = (B'*Ainv*B + C) \ (B'*Ainv*F + G);
    % q = Ainv*(F - B*u);
    Z = B'*Ainv*B + C;
    temp = Ainv*F;
    u = Z \ (B'*temp + G);
    q = temp - Ainv*(B*u);
    if show_spy
        spy(Z);
        print(gcf, '-depsc2', '-vector', 'spySS', '-r300');
    end
end

end