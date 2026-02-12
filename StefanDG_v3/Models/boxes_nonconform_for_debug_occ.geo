// Gmsh project created on Sat Jul 05 00:15:51 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 2, 1, 1};
//+
Point(10) = {1, 0, 0, 1.0};
//+
Point(11) = {1, 1, 0, 1.0};
//+
Point(12) = {1, 1, 1, 1.0};
//+
Point(13) = {1, 0, 1, 1.0};
//+
Line(13) = {12, 13};
//+
Line(14) = {13, 10};
//+
Line(15) = {10, 11};
//+
Line(16) = {11, 12};
//+
Curve Loop(7) = {13, 14, 15, 16};
//+
Plane Surface(7) = {7};
//+
//Surface{7} In Volume{1};
//+
MeshSize{ PointsOf{ Volume{1}; } } = 1;
//+
Physical Volume("Water", 17) = {1};
//+
Physical Surface("D", 18) = {1, 2};
