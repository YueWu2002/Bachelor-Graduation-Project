function [u, du, flux] = pLaplace1D_LDG(f, xx, gl, gr, p, p_equ, alpha, C12, type, plot_con)
% 1D p-Laplace equation with purely Dirichlet BC
% p: polynomial order
% p_equ: p-value
% alpha: C_{11} penalty order, 1 or 0


[qx, qw] = quadrature_nodal1D(p+1, 'legendre');     % quadrature
[nx, ~] = quadrature_nodal1D(p+1, type);            % collocation
[M, Minv, D, ~, ~, Tl, Tr] = matrix_nodal1D(nx);    % collocation
Pnq = interpolatemat(nx, qx);                       % maps nodal values at collocation nodes to those at quadrature nodes
Pqn = Minv*Pnq'*diag(qw);                           % quadrature-based interpolation operator (maps values at quadrature nodes to nodal representations)

q_equ = p_equ / (p_equ - 1);

N = length(xx) - 1; % number of elems
hh = xx(2:end) - xx(1:end-1);

Pnq_big = kron(speye(N), Pnq);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avg = [0.5*Tr, 0.5*Tl];
jump = [Tr, -Tl];
ETA = 1.0 * ones(N+1, 1);


% assemble QWEI (quadrature weights vector, based on quadrature points)
QWEI = zeros(N*length(qw), 1);
for k = 1: N
    ind = (k-1)*length(qw)+1 : k*length(qw);
    QWEI(ind) = (hh(k)/2) * qw(:);
end


% assemble JMP (maps nodal values to jumps, without Neumann boundaries)
JMP = sparse(N+1, N*(p+1));
for k = 1: N-1
    % interior jump between I_{k} and I_{k+1}
    ind = (k-1)*(p+1)+1 : (k+1)*(p+1);
    JMP(k+1, ind) = jump;
end
% add Dirichlet boundary jumps
JMP(1, 1:p+1) = -Tl;
JMP(N+1, (N-1)*(p+1)+1 : N*(p+1)) = Tr;


% assemble H_perp (characteristic height, without Neumann boundaries)
H_perp = zeros(N+1, 1);
for k = 1: N-1
    H_perp(k+1) = sqrt(hh(k) * hh(k+1));
end
H_perp(1) = hh(1);
H_perp(N+1) = hh(N);


% assemble Ainv (p=2) (inverse of the mass matrix)
Ainv_block = cell(N,1);
for k = 1: N
    Ainv_block{k} = sparse((2/hh(k)) * Minv);
end
Ainv = blkdiag(Ainv_block{1:N});


% assemble B
B_block = cell(N,1);
for k = 1: N
    B_block{k} = sparse(-((hh(k)/2) *M) * ((2/hh(k)) * D));
end
B = blkdiag(B_block{1:N});
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
    C11 = ETA(i+1) / H_perp(i+1)^alpha;
    C(ind, ind) = C(ind, ind) + jump'*C11*jump;
end
% add Dirichlet boundary penalty terms
C11 = ETA(1) / H_perp(1)^alpha;
ind = 1: p+1;
C(ind, ind) = C(ind, ind) + Tl'*C11*Tl;
C11 = ETA(N+1) / H_perp(N+1)^alpha;
ind = (N-1)*(p+1)+1 : N*(p+1);
C(ind, ind) = C(ind, ind) + Tr'*C11*Tr;


% assemble F (only contain Dirichlet's data)
F = sparse(N*(p+1), 1);
ind = 1: p+1;
F(ind) = F(ind) - gl*Tl';
ind = (N-1)*(p+1)+1 : N*(p+1);
F(ind) = F(ind) + gr*Tr';


% G = F_dual + D_pnlt_dual + N_dual

% assemble F_dual (the source as a dual fun)
F_dual = zeros(N*(p+1), 1);
for k = 1: N
    base = (k-1)*(p+1);
    xc = 0.5*(xx(k) + xx(k+1));

    % quadrature-based interpolation
    fvals = zeros(length(qx),1);
    for j = 1: length(qx)
        xp = xc + (hh(k)/2)*qx(j);
        fvals(j) = f(xp);
    end
    F_dual(base+1 : base+(p+1)) = Pnq'*((hh(k)/2)*diag(qw))*fvals;
end


% assemble D_pnlt_dual (Dirichlet penalty as a dual fun)
D_pnlt_dual = sparse(N*(p+1), 1);
C11 = ETA(1) / H_perp(1)^alpha;
ind = 1: p+1;
D_pnlt_dual(ind) = D_pnlt_dual(ind) + C11*gl*Tl';
C11 = ETA(N+1) / H_perp(N+1)^alpha;
ind = (N-1)*(p+1)+1 : N*(p+1);
D_pnlt_dual(ind) = D_pnlt_dual(ind) + C11*gr*Tr';


% assemble N_dual (Neumann as a dual fun)
N_dual = sparse(N*(p+1), 1);


% assemble the DG distributional gradient operator
GRAD_in = -Ainv*B;  % gradient operator with zero Dirichlet data
GRAD_dbc = Ainv*F;  % gradient vector with only Dirichlet data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u = (B'*Ainv*B + C) \ (B'*GRAD_dbc + F_dual + D_pnlt_dual + N_dual);   % starts with Poisson solution
E1 = compute_energy(N, QWEI, Pnq_big, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u, p_equ);
rho = 1.0;  % step length, should start with 1.0 (because it's best for the Poisson equation)
maxit = 10000;
it = 0;
should_continue = true;

if plot_con
    E = [E1]; % starts from zero
    rhos = [];
    grad_norms = [];
end

while should_continue && it < maxit
    [grad, grad_norm] = solve_grad(N, QWEI, Pnq_big, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u, p_equ, 2, 1e-16);

    fun = @(t) compute_energy(N, QWEI, Pnq_big, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u - t.*grad, p_equ);

    [rho, E1] = golden_search_oneside(fun, E1, rho, 1e-6 * rho); % 1e-6 relative error should be good
    u = u - rho.*grad;

    if plot_con
        E = [E, E1];
        rhos = [rhos, rho];
        grad_norms = [grad_norms, grad_norm];
    end

    if (abs(rho) < 1e-16) || (abs(rho) > 1e+8) || (grad_norm < 1e-16)
        should_continue = false;
        % in most cases, we stop here due to rho == 0 (this means that f'(0) >= 0. which might occur due to finite accuracy of floating point operations)
    end

    it = it + 1;
end

if plot_con
    figure(1);
    plot(0:length(E)-1, E);
    xlabel('iteration');
    ylabel('J_{h}[u_{h}]');
    grid on;

    figure(2);
    semilogy(rhos);
    xlabel('iteration');
    ylabel('step size');
    grid on;

    figure(3);
    semilogy(grad_norms);
    xlabel('iteration');
    ylabel('||gradient||');
    grid on;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post-processing

du = GRAD_dbc + GRAD_in * u; % can be modified to improve accuracy
du = Ainv*(F - B*u);
flux = zeros(N*(p+1), 1);
for k = 1: N
    ind = (k-1)*(p+1)+1 : k*(p+1);

    % quadrature-based interpolation
    flux(ind) = Pqn*Amap(Pnq*du(ind), p_equ);
end

end


function [E] = compute_energy(N, QWEI, Pnq_big, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u, p_equ)
% compute the energy corresponding to the p-Laplace LDG minimization form
% u must be a column-vector

du = GRAD_dbc + GRAD_in * u; % could be improved to enhance floating point stability
jmp_u = JMP*u + sparse([1,N+1], 1, [gl, -gr], N+1, 1);
E = (sum(QWEI .* abs(Pnq_big*du).^p_equ) ...
    + sum(ETA .* (H_perp.^alpha) .* abs(jmp_u./(H_perp.^alpha)).^p_equ)) / p_equ ...
    - dot(F_dual + N_dual, u);
end


function [grad, grad_norm] = solve_grad(N, QWEI, Pnq_big, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u, p_equ, type, delta)
% solve the steepest gradient at u, and its related norm

% D_DG[u, g_D] == Ainv * (F - B*u)
% GRAD_in = -Ainv*B;  % gradient operator with zero Dirichlet data
% GRAD_dbc = Ainv*F;  % gradient vector with only Dirichlet data
du = GRAD_dbc + GRAD_in * u;
jmp_u = JMP*u + sparse([1,N+1], 1, [gl,-gr], N+1, 1);

GRAD_in_nq = Pnq_big * GRAD_in;
du_q = Pnq_big * du;

switch type
case 1
    % under H01 norm
    MAT = GRAD_in_nq' * (QWEI .* GRAD_in_nq) ...
    + JMP' * ( (ETA./(H_perp.^alpha)) .* JMP);
case 2
    % under linearized W1p norm / weighted H01 norm
    MAT = (QWEI .* GRAD_in_nq)' * (correct_wei(du_q, p_equ, delta) .* GRAD_in_nq) ...
    + JMP' * ( ( (ETA./(H_perp.^alpha)) .* correct_wei(jmp_u./(H_perp.^alpha), p_equ, delta) ) .* JMP);
case 3 % experimental
    % under the local metric given by the functional J
    MAT = GRAD_in_nq' * (QWEI .* (correct_wei(du_q, p_equ, delta) .* GRAD_in_nq)) ...
    + JMP' * ( ( (ETA./(H_perp.^alpha)) .* correct_wei(jmp_u./(H_perp.^alpha), p_equ, delta) ) .* JMP) ...
    + (p_equ - 2) * (...
        (QWEI .* (du_q .* GRAD_in_nq))' * (correct_fullwei(du_q, p_equ, delta) .* (du_q .* GRAD_in_nq)) ... 
        + (jmp_u .* JMP)' * ((ETA .* correct_fullwei(jmp_u./(H_perp.^alpha), p_equ, delta) .* (H_perp.^(-3*alpha)) .* jmp_u) .* JMP) ...
    );
case 4 % extremely slow
    % under standard L2 norm
    MAT = Pnq_big' * (QWEI .* Pnq_big);
otherwise
    error('Unknown preconditioner type!');
end

RHS = GRAD_in_nq' * (QWEI .* Amap(du_q, p_equ)) ...
+ JMP' * (ETA .* Amap(jmp_u./(H_perp.^alpha), p_equ)) ...
- F_dual - N_dual;

grad = MAT \ RHS;

grad_grad = GRAD_in * grad;
jmp_grad = JMP*grad;

grad_grad_q = Pnq_big * grad_grad;

switch type
case 1
    grad_norm = sqrt(dot(grad_grad_q, QWEI .* grad_grad_q) ...
    + dot(jmp_grad, (ETA ./ (H_perp.^alpha)) .* jmp_grad));
case 2
    grad_norm = sqrt( sum(correct_wei(du_q, p_equ, delta) .* QWEI .* abs(grad_grad_q).^2) + sum(correct_wei(jmp_u ./ (H_perp.^alpha), p_equ, delta) .* (ETA ./ (H_perp.^alpha)) .* abs(jmp_grad).^2));
case 3 % experimental
    grad_norm = sqrt( ...
        sum(correct_fullwei(du_q, p_equ, delta) .* QWEI .* (p_equ-1) .* abs(du_q .* grad_grad_q).^2) ...
        + sum(ETA .* correct_fullwei(jmp_u ./ (H_perp.^alpha), p_equ, delta) .* (p_equ-1) .* (H_perp.^(-3*alpha)) .* abs(jmp_u .* jmp_grad).^2) ...
    );
case 4 % extremely slow
    grad_q = Pnq_big * grad;
    grad_norm = sqrt( sum(QWEI .* abs(grad_q).^2) );
otherwise
    error('Unknown preconditioner type!');
end

end