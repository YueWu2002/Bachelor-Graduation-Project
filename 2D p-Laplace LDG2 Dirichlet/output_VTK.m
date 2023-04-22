function [] = output_VTK(mesh, sol, p, name)
% output the format in VTK format
% 
% intput: 
%   mesh:   mesh data
%   sol:    numerical solution
%   p:      polynomial order, p>=0
%   name:   output file name
%   title:  title in file
% 
% output: 
%   
% mesh should have at least the following structure: 
%            nodes: [3697×2 double]
%            elems: [7168×3 int32]
%            faces: [10864×2 int32]
%        elem2face: [7168×3 int32]
%     elem2facedir: [7168×3 logical]
%        face2elem: [10864×2 int32]
%         bd_faces: [224×1 int32]
%       face_types: [10864×1 int32]
%        elem2elem: [7168×3 int32]
% 
% sol should have at least the following structure: 
%        q: [110×7168 double]
%        u: [55×7168 double]
%        x: [55×7168 double]
%        y: [55×7168 double]

% open a file and append the content
fid = fopen(name, 'a');

% version information
fprintf(fid, '# vtk DataFile Version 3.0\n');

% title
fprintf(fid, 'vtk Output\n');

% storage format
fprintf(fid, 'ASCII\n');

% grid format (unstructured)
fprintf(fid, 'DATASET UNSTRUCTURED_GRID\n');
fprintf(fid, '\n');

% points data
fprintf(fid, 'POINTS %d double\n', size(mesh.nodes,1));
for k = 1: size(mesh.nodes,1)
    fprintf(fid, '%f %f 0.0\n', mesh.nodes(k,1:2));
end
fprintf(fid, '\n');

% elements' nodes
fprintf(fid, 'CELLS %d %d\n', [size(mesh.elems,1), (1+3)*size(mesh.elems,1)]);
for k = 1: size(mesh.elems,1)
    % warning: index starts from 0
    fprintf(fid, '3 %d %d %d\n', mesh.elems(k,1:3)-1);
end
fprintf(fid, '\n');

% elements' types
fprintf(fid, 'CELL_TYPES %d\n', size(mesh.elems,1));
for k = 1: size(mesh.elems,1)
    fprintf(fid, '5\n'); % 5 stands for triangle
end
fprintf(fid, '\n');


fclose(fid);

end