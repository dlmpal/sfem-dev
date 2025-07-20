//+
SetFactory("OpenCASCADE");
Rectangle(1) = {0, 0, 0, 1, 1, 0};
//+
Physical Curve("Left", 5) = {4};
//+
Physical Curve("Right", 6) = {2};
//+
Physical Curve("Bottom", 7) = {1};
//+
Physical Curve("Top", 8) = {3};
//+
Physical Surface("Solid", 9) = {1};
//+//+
Transfinite Curve {3, 1} = 10 Using Progression 1;
//+
Transfinite Curve {2, 4} = 10 Using Progression 1;
//+
Transfinite Surface {1} = {1, 2, 3, 4};
