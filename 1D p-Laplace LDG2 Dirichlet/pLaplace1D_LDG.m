function [q, u, du] = pLaplace1D_LDG(f, xx, gl, gr, p, p_equ, alpha, C12, type, plot_con)
% 1D p-Laplace equation with purely Dirichlet BC
% p: polynomial order
% p_equ: p-value
% alpha: C_{11} penalty order, 1 or 0

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

avg = [0.5*Tr, 0.5*Tl];
jump = [Tr, -Tl];

% assemble QWEI (quadrature weights vector, based on quadrature points)
QWEI = zeros(N*(p+1), 1);
for k = 1: N
    ind = (k-1)*(p+1)+1 : k*(p+1);
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

% assemble ETA (penalty coefficients at faces, without Neumann boundaries)
ETA = 1.0 * ones(N+1, 1);

% assemble H_perp (characteristic height, without Neumann boundaries)
H_perp = zeros(N+1, 1);
for k = 1: N-1
    H_perp(k+1) = sqrt(hh(k)* hh(k+1));
end
H_perp(1) = hh(1);
H_perp(N+1) = hh(N);

% assemble Ainv (p=2)
Ainv_block = cell(N,1);
for k = 1: N
    Ainv_block{k} = sparse((2/hh(k)) * sparse(Minv));
end
Ainv = blkdiag(Ainv_block{:});


% assemble B
S = cell(N,1);
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
C = JMP' * ( (ETA./(H_perp.^alpha)) .* JMP);


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
    for j = 1: p+1
        xp = xc + 0.5*hh(k)*qx(j);
        F_dual(base + j) = f(xp)*qw(j)*(hh(k)/2);
    end
end

% assemble D_pnlt_dual (Dirichlet penalty as a dual fun)
D_pnlt_dual = sparse(N*(p+1), 1);
C11 = ETA(1) * H_perp(1)^(-alpha);
ind = 1: p+1;
D_pnlt_dual(ind) = D_pnlt_dual(ind) + C11*gl*Tl';
C11 = ETA(N+1) * H_perp(N+1)^(-alpha);
ind = (N-1)*(p+1)+1 : N*(p+1);
D_pnlt_dual(ind) = D_pnlt_dual(ind) + C11*gr*Tr';

% assemble N_dual (Neumann as a dual fun)
N_dual = sparse(N*(p+1), 1);


% assemble the DG distributional gradient operator
GRAD_in = -Ainv*B;  % gradient operator with zero Dirichlet data
GRAD_dbc = Ainv*F;  % gradient vector with only Dirichlet data


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u1 = (B'*Ainv*B + C) \ (B'*GRAD_dbc + F_dual + D_pnlt_dual + N_dual);   % starts with Poisson solution
E = [compute_energy(N, QWEI, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u1, p_equ)];
rho = 1.0;  % step length, should start with 1.0 (because it's best for the Poisson equation)
rhos = [];
grad_norms = [];
list = ((sqrt(5)+1)/2).^linspace(-1,1,7); % coefficient for rho (should include 1.0) (golden section rule)
trial = cell(length(list));
Etrial = zeros(length(list),1);
should_continue = true;
maxit = 10000;
it = 0;
while should_continue
    u0 = u1;
    E0 = E(end);

    [grad, grad_norm] = solve_grad(N, p, QWEI, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u0, p_equ, 2);
    grad_norms = [grad_norms, grad_norm];

    if grad_norm < 1e-16
        break;
    end

    % line search
    flag = true;
    while flag
        if rho < 1e-16
            flag = false;
            should_continue = false;
        elseif rho > 1e+8
            flag = false;
            should_continue = false;
            disp('fail to converge!');
        end

        % try line search
        for k = 1: length(list)
            trial{k} = u0 - (list(k)*rho)*grad;
            Etrial(k) = compute_energy(N, QWEI, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, trial{k}, p_equ);
        end

        [mEtrial, id] = min(Etrial(:));
        if (mEtrial > E0)
            rho = 0.25*rho;
        else
            flag = false;
            u1 = trial{id(1)};
            E1 = Etrial(id(1));
            rho = list(id(1))*rho;
            E = [E, E1];
            rhos = [rhos, rho];
        end
    end

    if (~should_continue || it >= maxit)
        break;
    end
    it = it + 1;
end

if plot_con
    figure(1);
    plot(E);
    xlabel('iteration');
    ylabel('J[u]');
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

u = u1;
du = GRAD_dbc + GRAD_in * u;
q = Amap(du, p_equ);
end

function [y] = correct_wei(x, p_equ, delta)
% delta + abs(x).^(p-2) or (delta + abs(x)).^(p-2)
if p_equ == 2
    y = ones(size(x));
elseif p_equ > 2
    y = delta + (abs(x).^(p_equ-2));
else
    y = (delta + abs(x)).^(p_equ-2);
end
end

function [y] = Amap(x, p_equ)
% y = x^(p-1)*sgn(x), if p~=1
y = sign(x).*(abs(x).^(p_equ-1));
end

function [E] = compute_energy(N, QWEI, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u, p_equ)
% compute the energy corresponding to the p-Laplace LDG minimization form

du = GRAD_dbc + GRAD_in * u;
jmp_u = JMP*u + sparse([1,N+1], 1, [gl,-gr], N+1, 1);
E = (dot(abs(du).^p_equ, QWEI(:)) ...
    + dot(ETA .* H_perp.^alpha, abs(jmp_u./(H_perp.^alpha)).^p_equ))/p_equ ...
    - dot(F_dual + N_dual, u);
end

function [grad, grad_norm] = solve_grad(N, p, QWEI, JMP, ETA, H_perp, GRAD_in, GRAD_dbc, F_dual, N_dual, alpha, gl, gr, u, p_equ, type)
% solve the steepest gradient at u, and its related norm

% D_DG[u, g_D] == Ainv * (F - B*u)
% GRAD_in = -Ainv*B;  % gradient operator with zero Dirichlet data
% GRAD_dbc = Ainv*F;  % gradient vector with only Dirichlet data
du = GRAD_dbc + GRAD_in * u;
jmp_u = JMP*u + sparse([1,N+1], 1, [gl,-gr], N+1, 1);
delta = 1e-4;

switch type
case 1
    % under H01 norm
    MAT = GRAD_in' * (QWEI .* GRAD_in) ...
    + JMP' * ( (ETA./(H_perp.^alpha)) .* JMP);
case 2
    % under linearized W1p norm / weighted H01 norm
    MAT = GRAD_in'*((correct_wei(du, p_equ, delta).*QWEI).*GRAD_in) ...
    + JMP' * ( ( (ETA./(H_perp.^alpha)) .* correct_wei(jmp_u./(H_perp.^alpha), p_equ, delta) ) .* JMP);
case 3
    % under standard L2 norm
    MAT = sparse(1:N*(p+1), 1:N*(p+1), QWEI, N*(p+1), N*(p+1));
end

RHS = GRAD_in'*(QWEI.*Amap(du, p_equ)) ...
+ JMP'*(ETA.*Amap(jmp_u./(H_perp.^alpha), p_equ)) ...
- F_dual - N_dual;

grad = MAT \ RHS;

grad_grad = GRAD_in * grad;
jmp_grad = JMP*grad;

switch type
case 1
    grad_norm = sqrt(dot(grad_grad, QWEI.*grad_grad) ...
    + dot(jmp_grad, (ETA./(H_perp.^alpha)).*jmp_grad));
case 2
    grad_norm = sqrt(dot(correct_wei(du, p_equ, delta).*grad_grad, QWEI.*grad_grad) ...
    + dot(correct_wei(jmp_u./(H_perp.^alpha), p_equ, delta).*jmp_grad, (ETA./(H_perp.^alpha)).*jmp_grad));
case 3
    grad_norm = sqrt(dot(grad, QWEI.*grad));
end
end