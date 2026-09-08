SetFactory("OpenCASCADE");

// =======================
// Parameters
// =======================
H     = 2.0;     // inlet height
h     = 1.0;     // outlet height
Lin   = 5.0;     // inlet length
Lout  = 20.0;    // outlet length

nx_in  = 40;
nx_out = 160;
ny_h   = 20;
ny_H   = 40;

structured = 1;   // 1 = structured quads, 0 = unstructured tris
lc = 0.4;         // mesh size for unstructured

// =======================
// Points
// =======================
// Inlet (left)
Point(1) = {0, 0, 0, lc};
Point(2) = {Lin, 0, 0, lc};
Point(3) = {Lin, -h, 0, lc};
Point(4) = {Lin+Lout, -h, 0, lc};
Point(5) = {Lin+Lout,  h, 0, lc};
Point(6) = {0, h, 0, lc};

// =======================
// Boundary lines
// =======================
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

// =======================
// Surface
// =======================
Curve Loop(1) = {6, 1, 2, 3, 4, 5};
Plane Surface(1) = {1};

// =======================
// Physical groups
// =======================
Physical Curve("Inlet", 7) = {6};
Physical Curve("Outlet", 8) = {4};
Physical Curve("Walls", 9) = {5, 1, 3, 2};
Physical Surface("Fluid", 10) = {1};
