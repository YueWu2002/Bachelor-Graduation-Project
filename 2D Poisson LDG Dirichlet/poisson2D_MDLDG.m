function [sol] = poisson2D_MDLDG(mesh, f, g_D, g_N, p, vel, show_spy, solve_method)
% [sol] = poisson2D_MDLDG(mesh, f, g_D, g_N, p, vel, show_spy)
%   2D Poisson equation with Dirichlet BC
% 
% input: 
%   mesh:   self-defined struct for a triangular mesh
%   f:      source term
%   g_D, g_N:   function handles for boundary conditions
%   p:      polynomial order
%   vel:    velocity for determining C_12 in the md-LDG formulation (constant non-zero vector for convenience)
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

% interior quadrature (for the ref element: the right triangle)
% PROBLEM: the assembly process takes too much time, we should parallelize it!
% WARNING: for better speed, please DO NOT use cells to construct block matrices! Cells are slow!
tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare matrices

Np = ((p+1)*(p+2))/2;

% edge matrices (for the ref element: [-1,1]) (checked.)
[qxe, ~] = gausslobatto(p+1);
[Me, ~, ~, ~, ~] = matrix_nodal1D(qxe);

% interior matrices (for the ref element: the right triangle) (checked.)
[xn, yn, zn] = node_lobatto2D(p, false);
[M, Minv, Dx, Dy, Pmn, ~] = matrix_nodal2D(xn, yn, zn, p);
[T] = trace_lobatto_nodal2D(p);

% interior quadrature (for the ref element: the right triangle)
[qx, qy, qz, qw, nquad] = quadrature_nodal2D(2*p, false); % algebraic exactness = 2*p
Pmq = zeros(nquad, Np); % maps Taylor basis coefficients to quadrature nodal values
for k = 1: nquad
    Pmq(k,:) = Taylor_basis_eval(qx(k), qy(k), p);
end
Pnq = Pmq / Pmn;  % maps nodal values to quadrature nodal values (for quadrature-based interpolation)
% Pnq = Pmq*Pnm;



% REMARK: we now do not use quadratures for volume integrals, instead we interpolate the function at nodal points and use the mass matrix for integration

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
face_nh = zeros(nface, 1);
face_len = zeros(nface, 1);
normals = zeros(nface, 2);
C11 = zeros(nface, 1);
C12 = zeros(nface, 2);
parfor k = 1: nface
    face_vec = mesh.nodes(mesh.faces(k,2), :) - mesh.nodes(mesh.faces(k,1), :);
    face_len(k) = norm(face_vec);
    nm_vec = [0, 1; -1, 0] * face_vec';
    normals(k,:) = nm_vec / norm(nm_vec);   % row vector
    
    % compute characteristic "h"
    if any(mesh.face2elem(k,:) == 0)
        % boundary face (use normal length)

        % face h
        face_nh(k) = 2*areas(max(mesh.face2elem(k,:))) / face_len(k);

        % Dirichlet BC, we should use positive C11
        C11(k) = 10.0 / face_nh(k);

        % C12 need not to be defined on boundaries
        C12(k,:) = 0;
    else
        % interior face

        % face h (use geometric averatge of normal length)
        face_nh(k) = sqrt( (2*areas(mesh.face2elem(k,1)) / face_len(k)) * (2*areas(mesh.face2elem(k,2)) / face_len(k)));

        % use the minimal-dissipation LDG penalty method
        C11(k) = 0;
        C12(k,:) = 0.5 * normals(k,:) * modified_sign(dot(vel, normals(k,:)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% assemble matrices (A(block-diagonal), B, C, F, G)
% define the variables: [q; u], length = Nelem*(2*Np) + Nelem*Np, where in "each element" q = [qx; qy]
% arrangement: [qx(1); qy(1); qx(2); qy(2); ...; qx(nelem); qy(nelem); u(1); ...; u(nelem)]

A = sparse(nelem*(2*Np), nelem*(2*Np));
Ainv = sparse(nelem*(2*Np), nelem*(2*Np));
B = sparse(nelem*(2*Np), nelem*Np);
C = sparse(nelem*Np, nelem*Np);
F = zeros(nelem*(2*Np), 1);
G = zeros(nelem*Np, 1);

% except for really-block-diagonal matrices, please do not use cells for assembly!
% blocks for volume integrals
A_block = cell(nelem,1);
Ainv_block = cell(nelem,1);
B_block = cell(nelem,1);
G_block = cell(nelem,1);

% functions for computing the node indices from an element index
ind_q = @(eid) (eid-1)*(2*Np)+1 : eid*(2*Np);   % row vector
ind_u = @(eid) (eid-1)*Np+1 : eid*Np;           % row vector


% (1) assemble volume integrals; (2) assemble face integrals (main time cost)

% element-loop for volume integrals (MUST be carried-out FIRST!!!)
parfor k = 1: nelem
    xx = mesh.nodes(mesh.elems(k,1:3), 1);
    yy = mesh.nodes(mesh.elems(k,1:3), 2);
    Jac = [
        xx(1) - xx(3), xx(2) - xx(3);
        yy(1) - yy(3), yy(2) - yy(3)
    ];
    M_elem = det(Jac)*M;
    Minv_elem = (det(Jac)^(-1))*Minv;

    % for A, block size = [2*Np, 2*Np]
    A_block{k} = sparse(blkdiag(M_elem, M_elem));
    Ainv_block{k} = sparse(blkdiag(Minv_elem, Minv_elem));

    % for B, block size = [2*Np, Np]
    B_block{k} = sparse(([Dx, Dy]*kron(inv(Jac), speye(Np)))' * M_elem);


    % for G, block size = [Np, 1]
    %{
    % pure interpolation
    fval = zeros(Np, 1);
    for i = 1: Np
        xp = dot(xx, [xn(i), yn(i), zn(i)]);
        yp = dot(yy, [xn(i), yn(i), zn(i)]);
        fval(i) = f(xp, yp);
    end
    G_block{k} = M_elem*fval;
    %}
    % quadrature-based interpolation
    fvals = zeros(nquad, 1);
    for i = 1: nquad
        xp = dot(xx, [qx(i), qy(i), qz(i)]);
        yp = dot(yy, [qx(i), qy(i), qz(i)]);
        fvals(i) = f(xp, yp);
    end
    G_block{k} = Pnq' * ((det(Jac)*qw) .* fvals);
end
A = blkdiag(A_block{1:nelem});
Ainv = blkdiag(Ainv_block{1:nelem});
B = blkdiag(B_block{1:nelem});
for k = 1: nelem
    G(ind_u(k)) = G_block{k};
end
% release temporary variables
clear A_block Ainv_block B_block G_block;


% REMARK: the face loop costs most of the time in the assembly process
% face loop for face integrals (MUST be carried-out after volume integrals)

% a method that is VERY time-efficient in the serial mode, however, the storage-cost is twice as needed
% this method could also be parallelized
iB = zeros((4*Np)*(2*Np), nface, 'int32');
jB = zeros((4*Np)*(2*Np), nface, 'int32');
vB = zeros((4*Np)*(2*Np), nface);
iC = zeros((2*Np)*(2*Np), nface, 'int32');
jC = zeros((2*Np)*(2*Np), nface, 'int32');
vC = zeros((2*Np)*(2*Np), nface);
iF = zeros(2*Np, nface, 'int32');
vF = zeros(2*Np, nface);
iG = zeros(Np, nface, 'int32');
vG = zeros(Np, nface);
parfor k = 1: nface
    [iB_new,jB_new,vB_new, iC_new,jC_new,vC_new, iF_new,vF_new, iG_new,vG_new] = face_integrals_assmebly(k, mesh, g_D, g_N, p, Np, qxe, Me, T, normals(k,:), C11(k), C12(k,:), face_len(k));

    iB(:,k) = fill_from_front(iB_new, iB(:,k));
    jB(:,k) = fill_from_front(jB_new, jB(:,k));
    vB(:,k) = fill_from_front(vB_new, vB(:,k));

    iC(:,k) = fill_from_front(iC_new, iC(:,k));
    jC(:,k) = fill_from_front(jC_new, jC(:,k));
    vC(:,k) = fill_from_front(vC_new, vC(:,k));

    iF(:,k) = fill_from_front(iF_new, iF(:,k));
    vF(:,k) = fill_from_front(vF_new, vF(:,k));

    iG(:,k) = fill_from_front(iG_new, iG(:,k));
    vG(:,k) = fill_from_front(vG_new, vG(:,k));
end
iB = iB(:); jB = jB(:); vB = vB(:); Bind = (iB ~= 0);
iC = iC(:); jC = jC(:); vC = vC(:); Cind = (iC ~= 0);
iF = iF(:); vF = vF(:); Find = (iF ~= 0);
iG = iG(:); vG = vG(:); Gind = (iG ~= 0);
B = B + sparse(iB(Bind), jB(Bind), vB(Bind), nelem*(2*Np), nelem*Np);
C = C + sparse(iC(Cind), jC(Cind), vC(Cind), nelem*Np, nelem*Np);
F = F + sparse(iF(Find), int32(1), vF(Find), nelem*(2*Np), 1);
G = G + sparse(iG(Gind), int32(1), vG(Gind), nelem*Np, 1);
clear iB jB vB Bind iC jC vC Cind iF jF Find iG jG Gind;

assemble_time = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve

tic;

switch solve_method
case {'direct', 'full'}
    ZZ = [A, B; B', -C];
    RHS = [F; -G];
    SOL = ZZ \ RHS;
    if show_spy
        spy(ZZ);
        print(gcf, '-depsc2', '-vector', 'spyZZ', '-r300');
    end
case 'reduce'
    % preferred method of solving
    SS = B'*Ainv*B + C;
    RHSS = B'*(Ainv*F) + G;
    u = SS \ RHSS;  % we have tested PCG, which is MUCH slower and MUCH less accurate!
    q = Ainv * (F - B*u);
    SOL = [q; u];
    if show_spy
        spy(SS);
        print(gcf, '-depsc2', '-vector', 'spySS', '-r300');
    end
otherwise
    disp('Error solve_method!');
end

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


function [iB,jB,vB, iC,jC,vC, iF,vF, iG,vG] = face_integrals_assmebly(k, mesh, g_D, g_N, p, Np, qxe, Me, T, normal, C11, C12, face_len)
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