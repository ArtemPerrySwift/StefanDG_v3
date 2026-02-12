// Gmsh project created on Sun Jan 25 00:15:00 2026
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
Physical Surface("C0", 25) = {7, 2};
//+
Translate {-1, 0, 0} {
  Volume{2}; 
}
//+
Box(3) = {3, 0, 0, 1, 1, 1};
//+
MeshSize{ PointsOf{ Volume{3}; } } = 0.5;
//+
Translate {-1, 0, 0} {
  Volume{3}; 
}
//+
Physical Volume("Water", 14) = {1, 2, 3};
//+
Physical Surface("D", 13) = {1, 14};
//+
Physical Surface("C1", 26) = {8, 13};
