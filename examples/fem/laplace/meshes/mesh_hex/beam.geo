SetFactory("OpenCASCADE");
Box(1) = {0, 0, 0, 5, 1, 1};
Box(2) = {5, 0, 0, 5, 1, 1};
Coherence;
Transfinite Curve{:} = 10;
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{:};
//+
Physical Surface("Left", 21) = {1};
//+
Physical Surface("Right", 22) = {7};
//+
Physical Volume("Solid", 23) = {1, 2};
//+
Physical Surface("Top", 24) = {4};
//+
Physical Surface("Bottom", 25) = {3};
