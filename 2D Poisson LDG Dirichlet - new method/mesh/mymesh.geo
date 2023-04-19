// Gmsh project created on Wed Mar 08 22:54:33 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {1, 1, 0, 2.0};
//+
Point(2) = {-1, 1, 0, 2.0};
//+
Point(3) = {-1, -1, 0, 2.0};
//+
Point(4) = {0, -1, 0, 2.0};
//+
Point(5) = {0, 0, 0, 2.0};
//+
Point(6) = {1, 0, 0, 2.0};
//+
Point(7) = {0, 1, 0, 1.0};
//+
Point(8) = {-1, 0, 0, 1.0};
//+
Line(1) = {1, 6};
//+
Line(2) = {6, 4};
//+
Line(3) = {4, 3};
//+
Line(4) = {3, 8};
//+
Line(5) = {8, 2};
//+
Line(6) = {2, 7};
//+
Line(7) = {7, 1};
//+
Line(8) = {7, 5};
//+
Line(9) = {5, 4};
//+
Line(10) = {5, 8};
//+
Line(11) = {5, 6};
//+
Curve Loop(1) = {6, 8, 10, 5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {7, 1, -11, -8};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {10, -4, -3, -9};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 2, -9};
//+
Plane Surface(4) = {4};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Transfinite Surface {3};
