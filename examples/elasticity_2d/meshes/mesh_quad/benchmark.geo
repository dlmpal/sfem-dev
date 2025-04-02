//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 0.5 * 1e-3, 10e-3, 0};
//+
Physical Surface("Solid", 5) = {1};
//+
Physical Curve("Fixed", 6) = {1};
//+
Physical Curve("Right", 7) = {2};
//+
Physical Curve("Left", 8) = {4};
//+
Physical Curve("Upper", 9) = {3};
