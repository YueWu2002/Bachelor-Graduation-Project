// Gmsh project created on Tue Apr 25 10:47:15 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {1, 1, 0, 1.0};
//+
Point(2) = {1, 2, 0, 1.0};
//+
Point(3) = {2, 2, 0, 1.0};
//+
Point(4) = {2, 1, 0, 1.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 2};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
MeshSize {1} = 0.4;
