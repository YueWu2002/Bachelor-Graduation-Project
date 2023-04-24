function output_vtk(mesh, sol, p, name, only_plot_mesh)
% output the result in vtk format 
%   (only valid if the nodal set is the Lobatto set)
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
% 
% sol should have at least the following structure: 
%        q: [110×7168 double]
%        u: [55×7168 double]
%        x: [55×7168 double]
%        y: [55×7168 double]

if p < 0.0
    error('otuout_vtk: error p value!');
end

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

if p == 0 || only_plot_mesh == true
    % use the finite-volume type method to visualize

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

    if only_plot_mesh == false
        % output in the finite-volume style

        fprintf(fid, 'CELL_DATA %d\n', size(mesh.elems,1));
        fprintf(fid, '\n');

        % output u
        fprintf(fid, 'SCALARS u double 1\n');
        fprintf(fid, 'LOOKUP_TABLE table1\n');
        for k = 1: size(mesh.elems,1)
            fprintf(fid, '%f\n', sol.u(1,k));
        end
        fprintf(fid, '\n');
    end

elseif p > 0
    % use refinement to visualize

    Np = ((p+1)*(p+2))/2;

    % points data (all the output nodes, containing element vertices since we are using Lobatto nodes)
    fprintf(fid, 'POINTS %d double\n', numel(sol.u));
    for k = 1: size(mesh.elems,1)
        for j = 1: Np
            fprintf(fid, '%f %f 0.0\n', [sol.x(j,k), sol.y(j,k)]);
        end
    end
    fprintf(fid, '\n');

    
    % elements' nodes (contained in those output nodes)
    fprintf(fid, 'CELLS %d %d\n', [size(mesh.elems, 1), (1+3) * size(mesh.elems, 1)]);
    for k = 1: size(mesh.elems,1)
        fprintf(fid, '3 %d %d %d\n', (k-1)*Np + [Np,p+1,1] - 1);
    end
    fprintf(fid, '\n');

    % elements' types
    fprintf(fid, 'CELL_TYPES %d\n', size(mesh.elems,1));
    for k = 1: size(mesh.elems,1)
        fprintf(fid, '5\n'); % 5 stands for triangle
    end
    fprintf(fid, '\n');
    

    if only_plot_mesh == false

        fprintf(fid, 'POINT_DATA %d\n', numel(sol.u));
        fprintf(fid, '\n');

        % output u
        fprintf(fid, 'SCALARS u double 1\n');
        fprintf(fid, 'LOOKUP_TABLE table1\n');
        for k = 1: size(mesh.elems,1)
            for j = 1: Np
                fprintf(fid, '%f\n', sol.u(j,k));
            end
        end
        fprintf(fid, '\n');

        % output q
        fprintf(fid, 'VECTORS q double\n');
        for k = 1: size(mesh.elems,1)
            for j = 1: Np
                fprintf(fid, '%f %f 0.0\n', [sol.q(j,1,k), sol.q(j,2,k)]);
            end
        end
        fprintf(fid, '\n');

    end
end

fclose(fid);

end