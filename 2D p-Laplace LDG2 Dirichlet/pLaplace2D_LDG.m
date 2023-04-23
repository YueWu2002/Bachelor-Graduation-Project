function [sol] = pLaplace2D_MDLDG(mesh, f, g_D, g_N, p, p_equ, vel, show_spy, basis_type)
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
areas = zeros(nelem, 1);
parfor k = 1: nelem
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
parfor k = 1: nface
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
        % C12(k,:) = 0.5 * normals(k,:) * modified_sign(dot(vel, normals(k,:)));
        C12(k,:) = 0.0;
    end

    C11(k) = ETA(k) / face_nh(k);
end
outward_normals = sparse(outward_normals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble matrices (A(block-diagonal), B, C, F, G)
% define the variables: [q; u], length = Nelem*(2*Np) + Nelem*Np, where in "each element" q = [qx; qy]
% arrangement: [qx(1); qy(1); qx(2); qy(2); ...; qx(nelem); qy(nelem); u(1); ...; u(nelem)]

% maps edge nodal values to edge quadrature nodal values
Pnqe_big = kron(speye(nface), Pnqe);

% maps volume nodal values to volume quadrature nodal values for Qh
PnqQh_big = kron(speye(nelem), kron(speye(2), Pnq));

% global quadrature weights on the elements for Qh
QWEI_Qh = kron(areas(:)./0.5, repmat(qw(:), [2,1]));

% global quadrature weights on the face
QWEI_face = kron(face_len(:)./2.0, qwe(:));


% scalar-to-scalar jump and average operators (maps scalar functions to scalar functions)
% we split the operators into interior ones and boundary ones (may depend on boundary conditions, so we just use the trace operator)
ii = zeros((p+1), Np, 3, nelem, 'int32');
jj = zeros((p+1), Np, 3, nelem, 'int32');
JumpInvv = zeros((p+1), Np, 3, nelem);
AvgInvv = zeros((p+1), Np, 3, nelem);
TraceDirichletvv = zeros((p+1), Np, 3, nelem);
TraceNeumannvv = zeros((p+1), Np, 3, nelem);
parfor k = 1: nelem
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
parfor k = 1: nface
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
parfor k = 1: nelem
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
SourceF = reshape(SourceF, [],1);
B = blkdiag(B_block{1:nelem});
clear B_block;

% mass matrix of edge funciton space
Me_big = kron((face_len(:)./2.0).* speye(nface), Me);

B = B - reshape( [reshape(reshape(repmat(normals(:,1)', [p+1,1]), [],1).*(JumpIn + TraceNeumann), [Np*nface*(p+1), nelem]); reshape(reshape(repmat(normals(:,2)', [p+1,1]),[],1).*(JumpIn + TraceNeumann), [Np*nface*(p+1), nelem])], [nface*(p+1), 2*Np*nelem] )' * Me_big * (AvgIn + reshape(repmat(dot(C12(:,1:2), normals(:,1:2), 2)', [p+1,1]), [],1) .* (JumpIn + TraceNeumann));

% scalar-to-scalar jump operator on both interior faces and Dirichlet faces (similar to the definition in the paper)
JumpScalar = JumpIn + TraceDirichlet; % checked.
JumpScalarQXE = Pnqe_big * JumpScalar;

C = JumpScalar' * (kron(C11(:), ones(p+1,1)) .* Me_big) * JumpScalar;

F = (Pnqe_big * reshape( [reshape(reshape(repmat(outward_normals(:,1)', [p+1,1]), [],1) .* TraceDirichlet, [Np*nface*(p+1), nelem]); reshape(reshape(repmat(outward_normals(:,2)', [p+1,1]), [],1) .* TraceDirichlet, [Np*nface*(p+1), nelem])], [nface*(p+1), 2*Np*nelem] ))' * (QWEI_face(:) .* DirichletQXE(:));

SourceNeumann = (Pnqe_big * TraceNeumann)' * (QWEI_face(:) .* NeumannQXE(:));
SourceDirichlet = (Pnqe_big * TraceDirichlet)' * (QWEI_face(:) .* (kron(C11(:), ones(p+1,1)).* DirichletQXE(:)));

clear Me_big Pnqe_big;

% inverse of mass matrix of Qh
invMassQh = kron((0.5./areas(:)).*speye(nelem), kron(speye(2), Minv));

% assemble the DG distributional gradient operator
GRAD_in = sparse(-invMassQh*B);  % gradient operator with zero Dirichlet data
GRAD_dbc = sparse(invMassQh*F);  % gradient vector with only Dirichlet data
GRAD_in_QXE = sparse(PnqQh_big * GRAD_in);
GRAD_dbc_QXE = sparse(PnqQh_big * GRAD_dbc);
clear invMassQh PnqQh_big;

assemble_time = toc;
% end assmebly process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin solve
tic;

% exact solve the Poisson equation
% u = (C - B'*GRAD_in) \ (B'*GRAD_dbc + (SourceF + SourceNeumann + SourceDirichlet));

% initial value: Poisson equation solved by one-step gradient descent
u = zeros(Np, nelem);
RES = compute_RESIDUE(p_equ, p, Np, nelem, nface, u, GRAD_in, GRAD_in_QXE, GRAD_dbc, JumpScalar, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, Pnqe, ETA, face_nh, qw, qwe, QWEI_Qh, QWEI_face);
u(:) = u(:) - (C - B'*GRAD_in) \ RES(:);

% to be continued: code the iterative part


q = GRAD_in*u(:) + GRAD_dbc;
solve_time = toc;
% end solve
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin post-processing

% nodal positions where solutions are provided
x = nan(Np, nelem);
y = nan(Np, nelem);
parfor k = 1: nelem
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    for i = 1: Np
        x(i, k) = dot(xx, [xn(i), yn(i), zn(i)]);
        y(i, k) = dot(yy, [xn(i), yn(i), zn(i)]);
    end
end

sol = struct(...
    'q', reshape(q, [2*Np, nelem]), ... 
    'flux', zeros(2*Np, nelem), ...
    'u', reshape(u, [Np, nelem]), ... 
    'x', x, ... 
    'y', y ... 
);

fprintf('assembly time = %f s, solving time = %f, NDOF = %d \n', assemble_time, solve_time, nelem*Np);

% end post-processing
end


function [RES] = compute_RESIDUE(p_equ, p, Np, nelem, nface, u, GRAD_in, GRAD_in_QXE, GRAD_dbc, JumpScalar, JumpScalarQXE, DirichletQXE, SourceF, SourceNeumann, Pnq, Pnqe, ETA, face_nh, qw, qwe, QWEI_Qh, QWEI_face)
% compute the residue in the p-Laplace equation (aka. Frechet derivative of the energy funcitonal J_h)
% checked only for p_equ == 2

% jump[u] at face collocation nodes, without information from Dirichlet BC
jump_u = reshape(full(JumpScalar * u(:)), [p+1, nface]);

% jump[u - g_D] at face quadrature nodes, with information from Dirichlet BC
jump_u_qxe = nan(length(qwe), nface);
parfor k = 1: nface
    jump_u_qxe(:,k) = Pnqe * jump_u(:,k) - DirichletQXE(:,k);
end

jump_u_penalty = nan(length(qwe), nface); % penalty terms arising from jump[u - g_D] at face quadrature nodes
parfor k = 1: nface
    temp = nan(length(qwe), 1);
    for i = 1: length(qwe)
        temp(i) = ETA(k) * A_map(jump_u_qxe(i,k) / face_nh(k), p_equ);
    end
    jump_u_penalty(:,k) = temp(:);
end

q = reshape(full(GRAD_in*u(:) + GRAD_dbc), [Np, 2, nelem]);
Aq_qxe = nan(length(qw), 2, nelem); % split the x and y components
parfor k = 1: nelem
    q_qxe_local = Pnq * q(:,:,k);
    Aq_qxe_local = nan(length(qw), 2);
    for i = 1: length(qw)
        Aq_qxe_local(i,:) = A_map(q_qxe_local(i,:), p_equ);
    end
    Aq_qxe(:,:,k) = Aq_qxe_local(:,:);
end

RES = GRAD_in_QXE' * (QWEI_Qh(:) .* Aq_qxe(:)) ... 
    + JumpScalarQXE' * (QWEI_face(:) .* jump_u_penalty(:)) ... 
    - SourceF(:) - SourceNeumann(:);

end