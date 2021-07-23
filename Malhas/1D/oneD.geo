// Gmsh project created on Sun Jun 11 13:28:27 2021
a = 1;
//+
Point(1) = {0, 0, 0, a};
//+
Point(2) = {10, 0, 0, a};
//+
Line(1) = {1, 2};
//+
//+
Physical Curve("line", 1) = {1};
//+
Physical Point("left", 3) = {1, 2};
//+
//Physical Point("right", 3) = {2};
