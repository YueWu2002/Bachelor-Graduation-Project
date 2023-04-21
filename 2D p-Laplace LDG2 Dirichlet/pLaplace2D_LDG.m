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

% PROBLEM: the assembly process takes too much time, we should parallelize it!
% WARNING: for better speed, please DO NOT use cells to construct block matrices! Cells are slow!

tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare matrices

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
% mesh preprocessing

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

        % face h
        face_nh(k) = 2*areas(eid) / face_len(k);

        if mesh.elem2facedir(eid,feid) == true
            outward_normals(k,:) = normals(k,:);
        else
            outward_normals(k,:) = -normals(k,:);
        end

        % for Dirichlet faces, we should use positive C11
        C11(k) = ETA(k) / face_nh(k);

        % for Neumann faces
        % ...(we have no Neumann faces)

        % C12 need not to be defined on boundaries
        C12(k,:) = 0.0;
    else
        % interior face

        % face h (use geometric averatge of normal length)
        face_nh(k) = min( 2*areas(mesh.face2elem(k,:)) / face_len(k) );

        C11(k) = ETA(k) / face_nh(k);
        % C12(k,:) = 0.5 * normals(k,:) * modified_sign(dot(vel, normals(k,:)));
        C12(k,:) = 0.0;
    end
end
outward_normals = sparse(outward_normals);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble matrices (A(block-diagonal), B, C, F, G)
% define the variables: [q; u], length = Nelem*(2*Np) + Nelem*Np, where in "each element" q = [qx; qy]
% arrangement: [qx(1); qy(1); qx(2); qy(2); ...; qx(nelem); qy(nelem); u(1); ...; u(nelem)]

% inverse of mass matrix of Qh
Ainv = kron((0.5./areas(:)).*speye(nelem), kron(speye(2), Minv));

% mass matrix of edge funciton space
Me_big = kron((face_len(:)./2) .* speye(nface), Me);

% maps edge nodal values to edge quadrature nodal values
Pnqe_big = kron(speye(nface), Pnqe);


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

        if mesh.face_types(faceid) > 0
            % boundary face

            % for Dirichlet faces
            if mesh.elem2facedir(k,j) == true
                TraceDirichletvv(:,:,j,k) = T(:,:,j);
            else
                TraceDirichletvv(:,:,j,k) = T(end:-1:1,:,j); % here, we don't use the signed version as in the jump operator
            end

            % for Neumann faces
            % ...(in this program, all BCs are Dirichlet)

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
        
        ii(:,:,j,k) = (faceid-1)*int32(p+1) + repmat((int32(1):int32(p+1))', [1, Np]);
        jj(:,:,j,k) = (k-1)*int32(Np) + repmat(int32(1):int32(Np), [p+1, 1]);

    end
end
JumpIn = sparse(ii(:), jj(:), JumpInvv(:), nface*(p+1), nelem*Np); % checked.
AvgIn = sparse(ii(:), jj(:), AvgInvv(:), nface*(p+1), nelem*Np); % seems to be OK
TraceDirichlet = sparse(ii(:), jj(:), TraceDirichletvv(:), nface*(p+1), nelem*Np); % checked.
TraceNeumann = sparse(ii(:), jj(:), TraceNeumannvv(:), nface*(p+1), nelem*Np); % checked.
clear ii jj JumpInvv AvgInvv TraceDirichletvv TraceNeumannvv;


% boundary condition as a function on the boundary (checked.)
SourceDirichlet = zeros(p+1, nface);
SourceNeumann = zeros(p+1, nface);
parfor k = 1: nface
    if mesh.face_types(k) > 0
        % boundary face

        xx = mesh.nodes(mesh.faces(k,1:2), 1);
        yy = mesh.nodes(mesh.faces(k,1:2), 2);
        
        % Dirichlet face
        % quadrature-based interpolation
        gval = zeros(length(qwe), 1);
        for i = 1: length(qwe)
            xp = dot([0.5*(1-qxe(i)), 0.5*(1+qxe(i))], xx);
            yp = dot([0.5*(1-qxe(i)), 0.5*(1+qxe(i))], yy);
            gval(i) = g_D(xp, yp);
        end
        SourceDirichlet(:,k) = Pqne * gval;

        % Neumann face
        % ... (no Neumann face)

    end
end
SourceDirichlet = sparse(reshape(SourceDirichlet, [], 1));
SourceNeumann = sparse(reshape(SourceNeumann, [], 1));


% assemble B and SourceF
B_block = cell(nelem,1);
SourceF = nan(Np,nelem);
parfor k = 1: nelem
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    detJac = det(Jac);
    invJac = inv(Jac);

    B_block{k} = sparse(([Dx, Dy] * kron(invJac, speye(Np)))' * (detJac.*M)); % checked.
    % to be optimized...

    % quadrature-based interpolation (checked.)
    fvals = zeros(nquad, 1);
    for i = 1: nquad
        xp = dot(xx, [qx(i), qy(i), qz(i)]);
        yp = dot(yy, [qx(i), qy(i), qz(i)]);
        fvals(i) = f(xp, yp);
    end
    SourceF(:,k) = Pnq' * ((det(Jac)*qw) .* fvals);

end
SourceF = reshape(SourceF, [], 1);
B = blkdiag(B_block{1:nelem});
clear B_block;
B = B - reshape( [reshape(kron(normals(:,1), ones(p+1,1)) .* JumpIn, [Np*nface*(p+1), nelem]); reshape(kron(normals(:,2), ones(p+1,1)) .* JumpIn, [Np*nface*(p+1), nelem])], [nface*(p+1), 2*Np*nelem] )' * Me_big * (AvgIn + kron(dot(C12(:,1:2), normals(:,1:2), 2), ones(p+1,1)) .* JumpIn);

% checked.
C = (JumpIn + TraceDirichlet)' * (kron(C11(:), ones(p+1,1)) .* Me_big) * (JumpIn + TraceDirichlet);

% checked.
F = reshape( [reshape(kron(outward_normals(:,1), ones(p+1,1)) .* TraceDirichlet, [Np*nface*(p+1), nelem]); reshape(kron(outward_normals(:,2), ones(p+1,1)) .* TraceDirichlet, [Np*nface*(p+1), nelem])], [nface*(p+1), 2*Np*nelem] )' * Me_big * SourceDirichlet;

% checked.
G = SourceF + TraceDirichlet' * (kron(C11(:), ones(p+1,1)) .* Me_big) * SourceDirichlet;


assemble_time = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve

tic;

% initial value: poisson
SS = B'*Ainv*B + C;
RHSS = B'*(Ainv*F) + G;
u = SS \ RHSS;  % we have tested PCG, which is MUCH slower and MUCH less accurate!









































q = Ainv * (F - B*u);
SOL = [q; u];

solve_time = toc;

sol = struct(...
    'q', SOL(1 : nelem*(2*Np)),...
    'u', SOL(nelem*(2*Np)+1 : nelem*(3*Np)),...
    'x', zeros(nelem*Np,1),...
    'y', zeros(nelem*Np,1)...
);
clear SOL;

for k = 1: nelem
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    for i = 1: Np
        sol.x((k-1)*Np+i) = dot(xx, [xn(i), yn(i), zn(i)]);
        sol.y((k-1)*Np+i) = dot(yy, [xn(i), yn(i), zn(i)]);
    end
end

fprintf('assembly time = %f s, solving time = %f, NDOF = %d \n', assemble_time, solve_time, nelem*(3*Np));

end


function [iB,jB,vB, iC,jC,vC, iF,vF, iG,vG] = face_integrals_assmebly(k, mesh, g_D, g_N, p, Np, qxe, Me, T, normal, C11, C12, face_len, Pqne)
% input: 
%   k:      index of the working face
%   others: shared variables
%   normal: unit normal vector for this face (in alignment with the face's orientation)
% 
% output:
%   triples [i,j,v] for assemble the B and C matrix (i,j,v are should be column vectors)
%   pairs [i,v] for assemble the F and G vector (i,v are column vectors)
% 
% REMARK: we proved by experiments that this function can greatly accelerate the face integral assembly process!

iB = []; jB = []; vB = [];
iC = []; jC = []; vC = [];
iF = []; vF = [];
iG = []; vG = [];

% mass matrix for the face
M_face = (face_len/2.0)*Me;

if all(mesh.face2elem(k,:) > 0)
    % interior face

    % face length for constructing the mass matrix on the edge
    % len = norm(mesh.nodes(mesh.faces(k,1), :) - mesh.nodes(mesh.faces(k,2), :));
    
    % element id
    eid1 = double(mesh.face2elem(k,1));
    eid2 = double(mesh.face2elem(k,2));

    % face id related to the element in [1,2,3]
    fide1 = find(mesh.elem2face(eid1,:) == k);
    fide2 = find(mesh.elem2face(eid2,:) == k);

    % face direction related to the element
    fdire1 = mesh.elem2facedir(eid1, fide1);
    fdire2 = mesh.elem2facedir(eid2, fide2);
    
    % outward unit normal vector, and sequence of boundary DOFs
    if fdire1 == true
        onm1 = normal;
        seq1 = 1 : 1 : p+1;
    else
        onm1 = -normal;
        seq1 = p+1 : -1 : 1;
    end
    if fdire2 == true
        onm2 = normal;
        seq2 = 1 : 1 : p+1;
    else
        onm2 = -normal;
        seq2 = p+1 : -1 : 1;
    end
    if fdire1 == fdire2
        disp("two neighbouring elements share a face with the same direction!");
    end

    % firstly, assemble the trace operator for scalar variables on this face
    trace1 = sparse(T(seq1, :, fide1));
    trace2 = sparse(T(seq2, :, fide2));

    % secondly, assemble the normal trace operator for vector variables on this face
    tracenm1 = kron(onm1, trace1);
    tracenm2 = kron(onm2, trace2);

    % for B: %∫Γ𝑜({𝑢ℎ} + 𝐂12 ⋅ ⟦𝑢ℎ⟧) ⋅ ⟦𝜏ℎ⟧ d𝑠
    block = -[tracenm1, tracenm2]' * M_face * [(0.5 + dot(C12, onm1))*trace1, (0.5 + dot(C12, onm2))*trace2];
    % B([ind_q(eid1), ind_q(eid2)], [ind_u(eid1), ind_u(eid2)]) = B([ind_q(eid1), ind_q(eid2)], [ind_u(eid1), ind_u(eid2)]) + block;
    rows = {1:2*Np, 2*Np+1:4*Np};
    cols = {1:Np, Np+1:2*Np};
    row_add = {(eid1-1)*(2*Np), (eid2-1)*(2*Np)};
    col_add = {(eid1-1)*Np, (eid2-1)*Np};
    for j = 1: 2
        for i = 1: 2
            [is, js, vs] = find(block(rows{i}, cols{j}));
            iB = [iB; reshape(is + row_add{i}, [], 1)];
            jB = [jB; reshape(js + col_add{j}, [], 1)];
            vB = [vB; reshape(vs, [], 1)];
        end
    end


    
    % for C:∫Γ𝑜 𝐶11 ⟦𝑢ℎ⟧ ⋅ ⟦𝑣ℎ⟧ d𝑠 (we simplify to the scalar operator)
    block = [trace1, -trace2]' * (C11*M_face) * [trace1, -trace2];
    % C([ind_u(eid1), ind_u(eid2)], [ind_u(eid1), ind_u(eid2)]) = C([ind_u(eid1), ind_u(eid2)], [ind_u(eid1), ind_u(eid2)]) + block;
    rows = {1:Np, Np+1:2*Np};
    cols = {1:Np, Np+1:2*Np};
    row_add = {(eid1-1)*Np, (eid2-1)*Np};
    col_add = {(eid1-1)*Np, (eid2-1)*Np};
    for j = 1: 2
        for i = 1: 2
            [is, js, vs] = find(block(rows{i}, cols{j}));
            iC = [iC; reshape(is + row_add{i}, [], 1)];
            jC = [jC; reshape(js + col_add{j}, [], 1)];
            vC = [vC; reshape(vs, [], 1)];
        end
    end

else
    % boundary face

    % nodes for the face
    xx = mesh.nodes(mesh.faces(k,1:2), 1);
    yy = mesh.nodes(mesh.faces(k,1:2), 2);

    % element id
    eid = double(max(mesh.face2elem(k,:)));

    % face id related to the element in [1,2,3]
    fide = find(mesh.elem2face(eid,:) == k);

    % face direction related to the element
    fdire = mesh.elem2facedir(eid, fide);
    
    % outward unit normal vector
    if fdire == true
        onm = normal;
        seq = 1 : 1 : p+1;
    else
        onm = -normal;
        seq = p+1 : -1 : 1;
    end

    % firstly, assemble the trace operator for scalar on this face
    trace = sparse(T(seq, :, fide));

    % secondly, assemble the normal trace operator for vector on this face
    tracenm = kron(onm, trace);

    % for Dirichlet BC
    block = trace' * (C11*M_face) * trace;
    % C(ind_u(eid), ind_u(eid)) = C(ind_u(eid), ind_u(eid)) + block;
    [iC, jC, vC] = find(block);
    iC = reshape(iC + (eid-1)*Np, [], 1);
    jC = reshape(jC + (eid-1)*Np, [], 1);
    vC = reshape(vC, [], 1);


    gval = zeros(p+1, 1);
    for i = 1: p+1
        xp = dot([0.5*(1-qxe(i)), 0.5*(1+qxe(i))], xx);
        yp = dot([0.5*(1-qxe(i)), 0.5*(1+qxe(i))], yy);
        gval(i) = g_D(xp, yp);
    end
    gval = Pqne * gval;

    block = tracenm' * M_face * gval;
    % F(ind_q(eid)) = F(ind_q(eid)) + tracenm' * M_face * gval;
    [iF, ~, vF] = find(block);
    iF = reshape(iF + (eid-1)*(2*Np), [], 1);
    vF = reshape(vF, [], 1);

    block = trace' * (C11*M_face) * gval;
    % G(ind_u(eid)) = G(ind_u(eid)) + trace' * (C11*M_face) * gval;
    [iG, ~, vG] = find(block);
    iG = reshape(iG + (eid-1)*Np, [], 1);
    vG = reshape(vG, [], 1);

    % for Neumann BC
    % no ---
    %
    
end

end

function [y] = fill_from_front(x,y)
y(1:numel(x)) = x(:);
end