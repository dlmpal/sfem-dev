//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 10, 1, 0};
//+
Physical Curve("Left", 5) = {4};
//+
Physical Curve("Right", 6) = {2};
//+
Physical Curve("Bottom", 7) = {1};
//+
Physical Curve("Top", 8) = {3};
//+
Physical Surface("Volume", 9) = {1};
