%  Matlab mesh
% mymesh2, Created by Gmsh
% ASCII
clear msh;
msh.nbNod = 8;
msh.POS = [
1 1 0;
-1 1 0;
-1 -1 0;
0 -1 0;
0 0 0;
1 0 0;
0 1 0;
-1 0 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 1 6 0
 6 4 0
 4 3 0
 3 8 0
 8 2 0
 2 7 0
 7 1 0
 7 5 0
 5 4 0
 5 8 0
 5 6 0
];
msh.TRIANGLES =[
 2 7 8 0
 8 7 5 0
 7 1 5 0
 5 1 6 0
 5 8 4 0
 4 8 3 0
 4 5 6 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
 5 0
 6 0
 7 0
 8 0
];
