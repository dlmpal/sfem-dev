//+
SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 20, 1, 0.5};
//+
Physical Surface("Fixed", 13) = {1};
//+
Physical Volume("Solid", 14) = {1};
//+
Physical Volume("Solid", 14) += {1};
//+
Physical Volume("Solid", 14) += {1};
//+
Physical Volume("Solid", 14) += {1};
//+
Physical Surface("Free", 15) = {2};
//+
Physical Surface("Upper", 16) = {4};
