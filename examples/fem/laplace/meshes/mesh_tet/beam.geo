//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 20, 1, 0.5};
//+
Physical Surface("Left", 13) = {1};
//+
Physical Volume("Solid", 14) = {1};
//+
Physical Volume("Solid", 14) += {1};
//+
Physical Volume("Solid", 14) += {1};
//+
Physical Volume("Solid", 14) += {1};
//+
Physical Surface("Right", 15) = {2};
//+
Physical Surface("Top", 16) = {4};
//+
Physical Surface("Bottom", 17) = {3};
