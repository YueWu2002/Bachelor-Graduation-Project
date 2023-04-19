function [mesh] = read_mesh(msh, plot_flag)
% read the planar triangular mesh from gmsh .m output

% correct inverted triangles
for k = 1: size(msh.TRIANGLES, 1)
    if (det([
        1, msh.POS(msh.TRIANGLES(k,1),1), msh.POS(msh.TRIANGLES(k,1),2);
        1, msh.POS(msh.TRIANGLES(k,2),1), msh.POS(msh.TRIANGLES(k,2),2);
        1, msh.POS(msh.TRIANGLES(k,3),1), msh.POS(msh.TRIANGLES(k,3),2)
    ]) < 0)
        msh.TRIANGLES(k,[1,2]) = msh.TRIANGLES(k,[2,1]);
    end
end

if plot_flag
    model = createpde;
    geometryFromMesh(model, msh.POS(:,1:2)', msh.TRIANGLES(:,1:3)', msh.TRIANGLES(:,4)+1);

    % plot the geometry
    pdegplot(model.Geometry);
    print(gcf, '-depsc2', '-vector', 'geom', '-r300');

    % plot the mesh
    pdemesh(model.Mesh);
    title([num2str(size(msh.TRIANGLES, 1)), ' cells, ', num2str(size(msh.POS, 1)), ' nodes']);
    print(gcf, '-depsc2', '-vector', 'mesh', '-r300');
end

% triangular mesh struct
mesh = struct(...
    'nodes', [],...
    'elems', [],...
    'faces', [],...
    'elem2face', [],...
    'elem2facedir', [],...
    'face2elem', [],...
    'bd_faces', []...
);
mesh.nodes = msh.POS(:,1:2);
mesh.elems = int32(msh.TRIANGLES(:,1:3));
[mesh.faces, mesh.elem2face, mesh.elem2elem, mesh.elem2facedir, mesh.face2elem, mesh.bd_faces] = make_face(mesh.elems);



end

function [face, elem2face, elem2elem, elem2facedir, face2elem, boundary] = make_face(t)
% [face, elem2face, elem2facedir, face2elem, boundary] = make_face(t)
% make faces from t-arrays (triangle)
% 
% t:         size = [nelem, 3], with node of each triangle in counter-clock-wise order
% face:      size = [nface, 2], with node indices of each face stored in increasing order
% elem2face: size = [nelem, 3], elem2face(k, 1) is the index of face elem(k, 2:3)
% elem2elem: size = [nelem, 3], elem2elem(k, 1) is the index of neighbour element sharing face elem2face(k, 1)
% elem2facedir: size = [nelem, 3], elem2facedir(k, :) indicats whether this face is the same direction as the saved one
% face2elem: size = [nface, 2], stroes the left-side element index and right-side one respectively (0 for empty)
% boundary:  size = [nbd, 1], stores the indices of boundary faces

nelem = size(t, 1);
after = [2; 3; 1];  % for cyclic looping the triangle

% make face through sparse matrix
node_max = max(t, [], 'all');
mat = sparse(t(:,1), t(:,2), 1, node_max, node_max, nelem) + sparse(t(:,2), t(:,3), 1, node_max, node_max, nelem) + sparse(t(:,3), t(:,1), 1, node_max, node_max, nelem);
mat = triu(mat + mat');
face_ind = int32(find(mat));
nface = length(face_ind);
face = zeros(nface, 2, 'int32');
face2elem = zeros(nface, 2, 'int32');
[face(:,1), face(:,2)] = ind2sub(size(mat), face_ind);
boundary = int32(find(mat(face_ind) == 1));

elem2face = zeros(nelem, 3, 'int32');
elem2facedir = true(nelem, 3);
for k = 1: nelem
    for j = 1: 3
        p1 = t(k, after(j));
        p2 = t(k, after(after(j)));
        if p1 < p2
            ii = sub2ind(size(mat), p1, p2);
            elem2face(k, j) = find(face_ind == ii);
            elem2facedir(k, j) = true;
            face2elem(elem2face(k, j), 1) = k;
        else
            ii = sub2ind(size(mat), p2, p1);
            elem2face(k, j) = find(face_ind == ii);
            elem2facedir(k, j) = false;
            face2elem(elem2face(k, j), 2) = k;
        end
    end
end

elem2elem = zeros(nelem, 3, 'int32');
for k = 1: nelem
    for i = 1: 3
        fid = elem2face(k,i);   % choose a face
        other_eidf = find(face2elem(fid,:)~=k);
        elem2elem(k,i) = face2elem(fid, other_eidf);
    end
end

end