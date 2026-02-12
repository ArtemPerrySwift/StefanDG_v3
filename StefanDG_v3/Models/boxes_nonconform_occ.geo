// Gmsh project created on Thu May 01 21:26:14 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
MeshSize{ PointsOf{ Volume{1}; } } = 0.5;
//+
Box(2) = {2, 0, 0, 1, 1, 1};
//+
MeshSize{ PointsOf{ Volume{2}; } } = 1;
//+
Physical Volume("Water", 14) = {1, 2};
//+
Physical Surface("D", 13) = {1, 8};
//+
Physical Surface("C0", 25) = {7, 2};
//+
Translate {-1, 0, 0} {
  Volume{2}; 
}
