function [sol] = pLaplace2D_MDLDG(mesh, f, g_D, g_N, p, p_equ, vel, show_spy, basis_type, plot_con)
% [sol] = pLaplace2D_MDLDG(mesh, f, g_D, g_N, p, p_equ, vel, show_spy, basis_type)
%   2D p-Laplace equation with Dirichlet BC
% 
% input: 
%   mesh:   self-defined struct for a triangular mesh
%   f:      source term
%   g_D, g_N:   function handles for boundary conditions
%   p:      polynomial order
%   vel:    velocity for determining C_12 in the md-LDG formulation (constant non-zero vector for convenience)
%   basis_type: type of polynomial basis
%   show_spy:   whether to plot the sparsity pattern of the matrix
% 
% output: 
%   sol:    solution struct: store (q, u) in each element, (x, y) for each data point
% (arrangement: [qx(1); qy(1); qx(2); qy(2); ...; qx(nelem); qy(nelem); u(1); ...; u(nelem)])
%
% REMARK: "mesh" data structure example:
%   mesh = 
%   struct with fields:
% 
%            nodes: [12×2 double]
%            elems: [14×3 int32]
%            faces: [25×2 int32]
%        elem2face: [14×3 int32]
%        elem2elem: [14×3 int32]
%     elem2facedir: [14×3 logical]
%        face2elem: [25×2 int32]
%         bd_faces: [8×1 int32]
%       face_types: [25x1 int32]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare matrices
tic;

Np = ((p+1)*(p+2))/2;

% edge matrices (for the ref element: [-1,1]) (checked.)
[qxe, qwe] = quadrature_nodal1D(p+1, 'legendre'); % quadrature points
[nxe, ~] = gausslobatto(p+1);
[Me, Minve, ~, ~, ~] = matrix_nodal1D(nxe);
Pnqe = interpmat_nodal1D(nxe, qxe);
Pqne = Minve*Pnqe'*diag(qwe);

% interior matrices (for the ref element: the right triangle) (checked.)
[xn, yn, zn] = node_lobatto2D(p, false);
[M, Minv, Dx, Dy, ~, ~] = matrix_nodal2D(xn, yn, zn, p, basis_type);
[T] = trace_lobatto_nodal2D(p);

% interior quadrature (for the ref element: the right triangle)
[qx, qy, qz, qw, nquad] = quadrature_nodal2D(2*p, false); % algebraic exactness = 2*p
Pnq = interpmat_nodal2D(xn, yn, p, qx, qy, basis_type); % maps nodal values to quadrature nodal values (for quadrature-based interpolation)
Pqn = Minv*Pnq'*diag(qw); % quadrature-based interpolation operator (maps values at quadrature nodes to nodal representations, which is not necessarily the inverse of Pnq)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesh pre-processing

nelem = size(mesh.elems, 1);
nface = size(mesh.faces, 1);

% compute element areas (checked.)
areas = nan(nelem, 1);
for k = 1: nelem
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    areas(k) = 0.5*det(Jac);
end

% compute the face length
% compute unit normal vectors for faces (point right to the face) (checked.)
% compute characteristic "h" in the normal direction for edges (for the penalty term) (checked.)
% compute the penalty coefficients (checked.)
eta = 10.0;
ETA = eta .* ones(nface, 1);
face_nh = zeros(nface, 1);
face_len = zeros(nface, 1);
normals = zeros(nface, 2); % outward unit normal vector related to faces
outward_normals = zeros(nface, 2); % outward unit normal vector related to the domain
C11 = zeros(nface, 1);
C12 = zeros(nface, 2);
for k = 1: nface
    face_vec = mesh.nodes(mesh.faces(k,2), :) - mesh.nodes(mesh.faces(k,1), :);
    face_len(k) = norm(face_vec);
    nm_vec = [0, 1; -1, 0] * face_vec';
    normals(k,:) = nm_vec / norm(nm_vec);   % row vector
    
    if mesh.face_types(k) > 0
        % boundary face

        eid = max(mesh.face2elem(k,:));
        feid = find(mesh.elem2face(eid,:)==k);

        if mesh.elem2facedir(eid,feid) == true
            outward_normals(k,:) = normals(k,:);
        else
            outward_normals(k,:) = -normals(k,:);
        end

        face_nh(k) = 2*areas(eid) / face_len(k);
        C12(k,:) = 0.0;
    else
        % interior face
        face_nh(k) = min( 2*areas(mesh.face2elem(k,:)) / face_len(k) );
        C12(k,:) = 0.5 * normals(k,:) * modified_sign(dot(vel, normals(k,:)));
        % C12(k,:) = 0.0;
    end

    C11(k) = ETA(k) / face_nh(k);
end
outward_normals = sparse(outward_normals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble matrices (A(block-diagonal), B, C, F, G)
% define the variables: [q; u], length = Nelem*(2*Np) + Nelem*Np, where in "each element" q = [qx; qy]
% arrangement: [qx(1); qy(1); qx(2); qy(2); ...; qx(nelem); qy(nelem); u(1); ...; u(nelem)]


% scalar-to-scalar jump and average operators (maps scalar functions to scalar functions)
% we split the operators into interior ones and boundary ones (may depend on boundary conditions, so we just use the trace operator)
ii = zeros((p+1), Np, 3, nelem, 'int32');
jj = zeros((p+1), Np, 3, nelem, 'int32');
JumpInvv = zeros((p+1), Np, 3, nelem);
AvgInvv = zeros((p+1), Np, 3, nelem);
TraceDirichletvv = zeros((p+1), Np, 3, nelem);
TraceNeumannvv = zeros((p+1), Np, 3, nelem);
for k = 1: nelem
    for j = 1: 3
        faceid = int32(mesh.elem2face(k,j));
        ii(:,:,j,k) = (faceid-1)*int32(p+1) + repmat((int32(1):int32(p+1))', [1, Np]);
        jj(:,:,j,k) = (k-1)*int32(Np) + repmat(int32(1):int32(Np), [p+1, 1]);
        if mesh.face_types(faceid) > 0
            % boundary face
            if mesh.elem2facedir(k,j) == true
                TraceDirichletvv(:,:,j,k) = T(:,:,j);
            else
                TraceDirichletvv(:,:,j,k) = T(end:-1:1,:,j); % here, we don't use the signed version as in the jump operator
            end
        else
            % interior face
            if mesh.elem2facedir(k,j) == true
                JumpInvv(:,:,j,k) = T(:,:,j);
                AvgInvv(:,:,j,k) = 0.5*T(:,:,j);
            else
                JumpInvv(:,:,j,k) = -T(end:-1:1,:,j); % use the signed version which helps to reduce cost
                AvgInvv(:,:,j,k) = 0.5*T(end:-1:1,:,j);
            end
        end
    end
end
JumpIn = sparse(ii(:), jj(:), JumpInvv(:), nface*(p+1), nelem*Np);
AvgIn = sparse(ii(:), jj(:), AvgInvv(:), nface*(p+1), nelem*Np);
TraceDirichlet = sparse(ii(:), jj(:), TraceDirichletvv(:), nface*(p+1), nelem*Np);
TraceNeumann = sparse(ii(:), jj(:), TraceNeumannvv(:), nface*(p+1), nelem*Np);
clear ii jj JumpInvv AvgInvv TraceDirichletvv TraceNeumannvv;


% boundary condition
DirichletQXE = zeros(length(qwe), nface); % store Dirichlet BC values at face quadrature nodes
NeumannQXE = zeros(length(qwe), nface); % store Neumann BC values at face quadrature nodes
for k = 1: nface
    if mesh.face_types(k) > 0
        % boundary face
        xx = mesh.nodes(mesh.faces(k,1:2), 1);
        yy = mesh.nodes(mesh.faces(k,1:2), 2);

        % Dirichlet boundary
        % quadrature-based interpolation
        gval = zeros(length(qwe), 1);
        for i = 1: length(qwe)
            xp = dot([0.5*(1-qxe(i)), 0.5*(1+qxe(i))], xx);
            yp = dot([0.5*(1-qxe(i)), 0.5*(1+qxe(i))], yy);
            gval(i) = g_D(xp, yp);
        end
        DirichletQXE(:,k) = gval(:);

        % Neumann boundary
        % ...(no)

    end
end

% assemble B and SourceF
B_block = cell(nelem, 1); % we assemble B using the divergence form
SourceF = nan(Np, nelem);
for k = 1: nelem
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    detJac = det(Jac);
    invJac = inv(Jac);

    % we assemble B using the divergence form
    B_block{k} = sparse([
        invJac(1,1).*Dx + invJac(2,1).*Dy, ...
        invJac(1,2).*Dx + invJac(2,2).*Dy
    ]' * (detJac.*M));

    % quadrature-based interpolation
    fvals = zeros(nquad, 1);
    for i = 1: nquad
        xp = dot(xx, [qx(i), qy(i), qz(i)]);
        yp = dot(yy, [qx(i), qy(i), qz(i)]);
        fvals(i) = f(xp, yp);
    end
    SourceF(:,k) = Pnq' * ((det(Jac)*qw) .* fvals);
end
B = blkdiag(B_block{1:nelem});
clear B_block;

% mass matrix of edge funciton space
Me_big = kron((face_len(:)./2.0).* speye(nface), Me);

% maps edge nodal values to edge quadrature nodal values
Pnqe_big = kron(speye(nface), Pnqe);

% maps volume nodal values to volume quadrature nodal values for Qh
PnqQh_big = kron(speye(nelem), kron(speye(2), Pnq));

% global quadrature weights on the elements for Vh
QWEI_Vh = kron(areas(:)./0.5, qw(:));

% global quadrature weights on the elements for Qh
QWEI_Qh = kron(areas(:)./0.5, repmat(qw(:), [2,1]));

% global quadrature weights on the face
QWEI_face = kron(face_len(:)./2.0, qwe(:));

B = B - reshape( [reshape(reshape(repmat(normals(:,1)', [p+1,1]), [],1).*(JumpIn + TraceNeumann), [Np*nface*(p+1), nelem]); reshape(reshape(repmat(normals(:,2)', [p+1,1]),[],1).*(JumpIn + TraceNeumann), [Np*nface*(p+1), nelem])], [nface*(p+1), 2*Np*nelem] )' * Me_big * (AvgIn + reshape(repmat(dot(C12(:,1:2), normals(:,1:2), 2)', [p+1,1]), [],1) .* (JumpIn + TraceNeumann));

% scalar-to-scalar jump operator on both interior faces and Dirichlet faces (similar to the definition in the paper)
JumpScalar = JumpIn + TraceDirichlet; % checked.
JumpScalarQXE = Pnqe_big * JumpScalar;

C = JumpScalar' * (kron(C11(:), ones(p+1,1)) .* Me_big) * JumpScalar;

F = reshape( [reshape(reshape(repmat(outward_normals(:,1)', [p+1,1]), [],1) .* TraceDirichlet, [Np*nface*(p+1), nelem]); reshape(reshape(repmat(outward_normals(:,2)', [p+1,1]), [],1) .* TraceDirichlet, [Np*nface*(p+1), nelem])], [nface*(p+1), 2*Np*nelem] )' * (Pnqe_big' * (QWEI_face(:) .* DirichletQXE(:)));

SourceNeumann = reshape(TraceNeumann' * (Pnqe_big' * (QWEI_face(:) .* NeumannQXE(:))), [Np, nelem]);
% SourceDirichlet = reshape(TraceDirichlet' * (Pnqe_big' * (QWEI_face(:) .* (kron(C11(:), ones(length(qwe),1)).* DirichletQXE(:)))), [Np, nelem]);

clear Pnqe_big TraceNeumann;


% inverse of mass matrix of Qh
invMassQh = kron((0.5./areas(:)).*speye(nelem), kron(speye(2), Minv));

% assemble the DG distributional gradient operator
GRAD_in = sparse(-invMassQh*B);  % gradient operator with zero Dirichlet data
clear B;
GRAD_dbc = reshape(sparse(invMassQh*F), [2*Np, nelem]);  % gradient vector with only Dirichlet data
clear invMassQh F;
GRAD_in_QXE = sparse(PnqQh_big * GRAD_in);

% mass matrix of Qh
MassQh = kron((areas(:)./0.5).*speye(nelem), kron(speye(2), M));

% store the decomposition
% (LU(auto) and LDL has the best accuracy, while Chol is faster)
PoissonMat_decom = decomposition(C + GRAD_in' * MassQh * GRAD_in, 'chol', 'lower', 'CheckCondition', false);
clear C;

assemble_time = toc;
% end assmebly process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin solve
% checked for p_equ == 2
tic;

% initial value: Poisson preconditioner
u = zeros(Np, nelem);
[w, ~] = compute_preconditioned_gradient(p_equ, p, Np, nelem, nface, u, GRAD_in, GRAD_in_QXE, GRAD_dbc, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, Pnqe, areas, ETA, face_len, face_nh, qw, qwe, QWEI_Qh, QWEI_face, PoissonMat_decom, 1, 1e-10);
u(:) = u(:) - 1.0 .* w(:);

E1 = compute_energy(p_equ, p, Np, nelem, nface, u, GRAD_in, GRAD_dbc, JumpScalar, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, PnqQh_big, Pnqe, ETA, face_nh, face_len, qw, qwe, areas, QWEI_Vh, QWEI_face);
rho = 1.0;  % step length, should start with 1.0 (because it's best for the Poisson equation)
maxit = 10000;
it = 0;
should_continue = true;

if plot_con
    E = [E1]; % starts from the energy of the initial guess
    rhos = [];
    grad_norms = [];
end

while should_continue && it < maxit

    % remark: big delta is good for stability, while small delta is good for accuracy!
    delta = 1e-14;

    [grad, grad_norm] = compute_preconditioned_gradient(p_equ, p, Np, nelem, nface, u, GRAD_in, GRAD_in_QXE, GRAD_dbc, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, Pnqe, areas, ETA, face_len, face_nh, qw, qwe, QWEI_Qh, QWEI_face, PoissonMat_decom, 2, delta);
    % warning: the Poisson preconditioner is extremely slow!
    % very recommended: use the weighted preconditioner!

    fun = @(t) compute_energy(p_equ, p, Np, nelem, nface, u(:) - t.*grad(:), GRAD_in, GRAD_dbc, JumpScalar, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, PnqQh_big, Pnqe, ETA, face_nh, face_len, qw, qwe, areas, QWEI_Vh, QWEI_face);

    [rho, E1] = golden_search_oneside(fun, E1, rho, 1e-6 * rho); % 1e-6 relative error should be good
    u(:) = u(:) - rho.*grad(:);

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


% post-process the solution
q = reshape(GRAD_in*u(:) + GRAD_dbc(:), [Np, 2, nelem]);
flux = nan(Np, 2, nelem);
parfor k = 1: nelem
    q_qx_local = Pnq * q(:,:,k);

    % quadrature-based interpolation
    flux_qx_local = nan(length(qw), 2);
    for i = 1: length(qw)
        flux_qx_local(i,:) = A_map(q_qx_local(i,:), p_equ);
    end
    flux(:,:,k) = Pqn * flux_qx_local;
end

solve_time = toc;
% end solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin packing the solution

% nodal positions where solutions are provided
x = nan(Np, nelem);
y = nan(Np, nelem);
parfor k = 1: nelem
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    for i = 1: Np
        x(i,k) = dot(xx, [xn(i), yn(i), zn(i)]);
        y(i,k) = dot(yy, [xn(i), yn(i), zn(i)]);
    end
end

sol = struct(...
    'q', reshape(q, [Np, 2, nelem]), ... 
    'flux', reshape(flux, [Np, 2, nelem]), ...
    'u', reshape(u, [Np, nelem]), ... 
    'x', reshape(x, [Np, nelem]), ... 
    'y', reshape(y, [Np, nelem]) ... 
);

fprintf('assembly time = %f s, solving time = %f, NDOF = %d \n', assemble_time, solve_time, nelem*Np);

% end packing
end


function [E] = compute_energy(p_equ, p, Np, nelem, nface, u, GRAD_in, GRAD_dbc, JumpScalar, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, PnqQh_big, Pnqe, ETA, face_nh, face_len, qw, qwe, areas, QWEI_Vh, QWEI_face)
% compute the energy functional Jh for the p-Laplace equation

% we find that: 'parfor' is very slow compared to vectorization

jump_u_qxe = reshape(full(JumpScalarQXE * u(:) - DirichletQXE(:)), [length(qwe), nface]);
q_qx = reshape(PnqQh_big * full(GRAD_in*u(:) + GRAD_dbc(:)), [length(qw), 2, nelem]);

%{
E = (... 
    sum(QWEI_Vh(:) .* reshape(sum(q_qx(:,:,:).^2, 2).^(p_equ/2.0), [],1)) ... 
    + sum(kron(ETA(:).*face_nh(:), ones(length(qwe),1)) .* QWEI_face(:) .* abs(jump_u_qxe(:)./kron(face_nh(:), ones(length(qwe),1))).^p_equ) ) / p_equ ... 
    - dot(SourceF(:) + SourceNeumann(:), u(:));
%}

E = (... 
    sum( (areas(:)./0.5)' .* (qw(:) .* reshape(sum(q_qx(:,:,:).^2, 2).^(p_equ/2.0), [length(qw), nelem])), 'all') ... 
    + sum( ((ETA(:).*face_nh(:)) .* (face_len(:)./2.0))' .* (qwe(:) .* abs(jump_u_qxe(:,:)./face_nh(:)').^p_equ), 'all') ) / p_equ ... 
    - dot(SourceF(:) + SourceNeumann(:), u(:));
    
end


function [w, w_norm] = compute_preconditioned_gradient(p_equ, p, Np, nelem, nface, u, GRAD_in, GRAD_in_QXE, GRAD_dbc, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, Pnqe, areas, ETA, face_len, face_nh, qw, qwe, QWEI_Qh, QWEI_face, PoissonMat, type, delta)
% compute the preconditioned gradient at u, and its related norm

% jump[u] at face quadrature nodes, with information from Dirichlet BC
jump_u_qxe = reshape(JumpScalarQXE * u(:) - DirichletQXE(:), [length(qwe), nface]);

jump_u_qxe_penalty = nan(length(qwe), nface); % penalty terms arising from jump[u - g_D] at face quadrature nodes
parfor k = 1: nface
    % jump[u - g_D] at face quadrature nodes, with information from Dirichlet BC
    temp = nan(length(qwe), 1);
    for i = 1: length(qwe)
        temp(i) = ETA(k) * A_map(jump_u_qxe(i,k) / face_nh(k), p_equ);
    end
    jump_u_qxe_penalty(:,k) = temp(:);
end

q_qx = reshape( Pnq * reshape(full(GRAD_in*u(:) + GRAD_dbc(:)), [Np, 2*nelem]), [length(qw), 2, nelem]);

Aq_qxe = nan(length(qw), 2, nelem); % split the x and y components
parfor k = 1: nelem
    Aq_qxe_local = nan(length(qw), 2);
    for i = 1: length(qw)
        Aq_qxe_local(i,:) = A_map(q_qx(i,:,k), p_equ);
    end
    Aq_qxe(:,:,k) = Aq_qxe_local(:,:);
end

%{
RES = GRAD_in_QXE' * (QWEI_Qh(:) .* Aq_qxe(:)) ... 
    + JumpScalarQXE' * (QWEI_face(:) .* jump_u_qxe_penalty(:)) ... 
    - SourceF(:) - SourceNeumann(:);
%}

RES = GRAD_in_QXE' * reshape( (areas(:)./0.5)' .* reshape( qw(:) .* Aq_qxe(:,:,:), [2*length(qw), nelem] ), [],1) + JumpScalarQXE' * reshape( (face_len(:)./2.0)' .* (qwe(:) .* jump_u_qxe_penalty(:,:)), [],1) - SourceF(:) - SourceNeumann(:);

% compute the preconditioned gradient
switch type
case 1 % Poisson preconditioner (avoid the use of quadrature-based interpolation operator Pqn, which help to maintain floating-point accuracy in the case of high order)
    w = PoissonMat \ RES(:); % the matrix can be stored because it never changes

    w_jump_qxe = reshape(JumpScalarQXE * w(:), [length(qwe), nface]);
    w_grad_qxe = reshape(GRAD_in_QXE * w(:), [length(qw), 2, nelem]);
    w_norm = sum((areas(:)./0.5)' .* (qw(:) .* reshape(sum(w_grad_qxe(:,:,:).^2, 2), [length(qw), nelem])), 'all') + sum(((ETA(:) ./ face_nh(:)) .* (face_len(:)./2.0))' .* (qwe(:) .* w_jump_qxe(:,:).^2), 'all');

case 2 % weighted preconditioner
    % error('weighted preconditioner not supported yet.');

    qx_weight = repmat((areas(:)./0.5)' .* (qw(:) .* correct_wei(sqrt(reshape(sum(q_qx(:,:,:).^2, 2), [length(qw), nelem])), p_equ, delta)), [2,1]); % size = [2*length(qw), nelem]

    qxe_weight = ((ETA(:) ./ face_nh(:)) .* (face_len(:)./2.0))' .* (qwe(:) .* correct_wei(jump_u_qxe(:,:)./face_nh(:)', p_equ, delta)); % size = [length(qwe), nface]

    MAT = GRAD_in_QXE' * (qx_weight(:) .* GRAD_in_QXE) ... 
        + JumpScalarQXE' * (qxe_weight(:) .* JumpScalarQXE);
    w = MAT \ RES(:);

    w_grad_qxe = reshape(GRAD_in_QXE * w(:), [2*length(qw), nelem]);
    w_jump_qxe = reshape(JumpScalarQXE * w(:), [length(qwe), nface]);
    w_norm = sum( qx_weight(:,:) .* w_grad_qxe(:,:).^2, 'all') + sum( qxe_weight(:,:) .* w_jump_qxe(:,:).^2, 'all');

case 3 % Newton-Raphson preconditioner
    error('Newton-Raphson preconditioner not supported yet.');

otherwise
    error('compute_preconditioned_gradient: unknown preconditioner type!');
end

end