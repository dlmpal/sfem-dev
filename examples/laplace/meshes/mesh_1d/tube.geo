//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Line(1) = {1, 2};
//+
Physical Point("Inlet", 2) = {1};
//+
Physical Point("Outlet", 3) = {2};
//+
Physical Curve("Tube", 4) = {1};
