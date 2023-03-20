function [q, u, du] = pLaplace1D_LDG(f, xx, gl, gr, p, p_equ, alpha, C12, type)
% 1D p-Laplace equation with purely Dirichlet BC
% p: polynomial order
% p_equ: p-value
% alpha: C_{11} penalty order, -1 or 0

switch type
case 'legendre' % should use legendre for integration due to stability requirements
    [qx, qw] = gausslegendre(p+1);
case 'lobatto'
    [qx, qw] = gausslobatto(p+1);
end
[M, Minv, D, ~, ~, Tl, Tr] = matrix_nodal1D(qx);
q_equ = p_equ / (p_equ - 1);

N = length(xx) - 1; % number of elems
hh = xx(2:end) - xx(1:end-1);


% 注意!!!跳跃算子请严格按照定义来实现，否则会出错!!!
avg = [0.5*Tr, 0.5*Tl];
jump = [Tr, -Tl];

% assemble B
for k = 1: N
    S{k} = sparse(-((hh(k)/2) * sparse(M)) * ((2/hh(k)) * sparse(D)));
end
B = blkdiag(S{:});
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
eta = 1.0;
% add jump terms
for i = 1: N-1
    % between I_{i} and I_{i+1}
    ind = (i-1)*(p+1)+1 : (i+1)*(p+1);
    C11 = eta * sqrt(hh(i) * hh(i+1))^alpha;
    C(ind, ind) = C(ind, ind) + jump'*C11*jump;
end
% add Dirichlet boundary terms
C11 = eta * hh(1)^alpha;
ind = 1: p+1;
C(ind, ind) = C(ind, ind) + Tl'*C11*Tl;
C11 = eta * hh(N)^alpha;
ind = (N-1)*(p+1)+1 : N*(p+1);
C(ind, ind) = C(ind, ind) + Tr'*C11*Tr;


% assemble F
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
    for j = 1: p+1
        xp = xc + 0.5*hh(k)*qx(j);
        G(base + j) = f(xp)*qw(j)*(hh(k)/2);
    end
end
C11 = eta * hh(1)^alpha;
ind = 1: p+1;
G(ind) = G(ind) + C11*gl*Tl';
C11 = eta * hh(N)^alpha;
ind = (N-1)*(p+1)+1 : N*(p+1);
G(ind) = G(ind) + C11*gr*Tr';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

qinds = zeros(N, p+1, 'int32');
for k = 1: N
    qinds(k,:) = (k-1)*(p+1)+1 : k*(p+1);
end

% assemble Ainv (p=2)
for k = 1: N
    Ainv_block{k} = sparse((2/hh(k)) * sparse(Minv));
end
Ainv = blkdiag(Ainv_block{:});

% assemble A (p=2)
for k = 1: N
    A_block{k} = sparse((hh(k)/2) * sparse(M));
end
A = blkdiag(A_block{:});
MASS = A;
MASSINV = Ainv;
ZZ = [A, B; B', -C];
RHS = [F; -G];

% initial value
x1 = [A, B; B', -C] \ RHS;

% Picard iteration
maxit = 50;
for it = 1: maxit
    x0 = x1;

    % partial weight
    % WEI = sparse(1:N*(p+1), 1:N*(p+1), correct(x0(1:N*(p+1)), q_equ, 1e-14), N*(p+1), N*(p+1), N*(p+1));
    % ZZ(1:N*(p+1), 1:N*(p+1)) = MASS*WEI; % A = MASS*WEI;
    WEI = correct(x0(1:N*(p+1)), q_equ, 1e-14);
    ZZ(1:N*(p+1), 1:N*(p+1)) = MASS .* WEI(:)';
    x1 = ZZ \ RHS;

    % Through experiment, the Schur-complement approach is significantly less stable!
    % Schur-complement algorithm for solving the local problem
    % inverse of partial weight
    % WEI = sparse(1:N*(p+1), 1:N*(p+1), 1./sqrt(correct(q0, q_equ, 1e-15)), N*(p+1), N*(p+1), N*(p+1));
    % Ainv = WEI'*MASSINV*WEI;
    % SS = B'*Ainv*B + C;
    % u1 = SS \ (B'*Ainv*F + G);
    % q1 = Ainv * (F - B*u0);
end

% Newton-Raphson iteration (faster)
maxit = 50;
for it = 1: maxit
    x0 = x1;
    RES = [
        F - MASS*Amap(x0(1:N*(p+1)), q_equ) - B*x0(N*(p+1)+1:2*N*(p+1));
        -G - B'*x0(1:N*(p+1)) + C*x0(N*(p+1)+1:2*N*(p+1))
    ];
    % ZZ(1:N*(p+1), 1:N*(p+1)) = MASS * sparse(1:N*(p+1), 1:N*(p+1), (q_equ-1)*Amap(x0(1:N*(p+1)),q_equ-1), N*(p+1), N*(p+1), N*(p+1));
    ZZ(1:N*(p+1), 1:N*(p+1)) = MASS * sparse(1:N*(p+1), 1:N*(p+1), (q_equ-1)*correct(x0(1:N*(p+1)), q_equ, 0.0), N*(p+1), N*(p+1), N*(p+1));
    x1 = x0 + ZZ \ RES;
end

q = x1(1 : N*(p+1));
u = x1(N*(p+1)+1 : end);
du = Amap(q, q_equ);
end

function [y] = correct(x, p_equ, delta)
% (delta + abs(x)).^(p-2)
% (delta + abs(x)).^(p-2) should >= delta
% (delta + abs(x)).^(2-p) should >= delta
y = x;
y(abs(y) < delta) = delta;
y(abs(y) > 1/delta) = 1/delta;
y(y==0) = eps;
y = abs(y).^(p_equ-2);
y(abs(y) < delta) = delta;
y(abs(y) > 1/delta) = 1/delta;
end

function [y] = Amap(x, p_equ)
% y = x^(p-1)*sgn(x), if p~=1
y = sign(x).*(abs(x).^(p_equ-1));
end