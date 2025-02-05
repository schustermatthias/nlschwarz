delta = 0.1;
cl__1 = 1;
lc = 1;
lc2 = 1;
lc3 = 1;
lc_fine = lc/2;
lc_point= lc/(15);


// Omega_1
Point(12) = {0.2, 0.75, 0, lc};
Point(13) = {0.8, 0.75, 0, lc};
Point(14) = {0.5, 0.85, 0, lc};
Point(15) = {0.5, 0.65, 0, lc};

// Omega_2
Point(17) = {0.2, 0.25, 0, lc};
Point(18) = {0.8, 0.25, 0, lc};
Point(19) = {0.5, 0.35, 0, lc};
Point(20) = {0.5, 0.15, 0, lc};


//------------------------------------------------------------------------------
// OMEGA_3
Point(1) = {0, 0, 0,  lc};
Point(2) = {0, 1, 0,  lc};
Point(3) = {1, 1, 0,  lc};
Point(4) = {1, 0, 0,  lc};

// Omega_I
Point(8) = {1. + delta, -delta, -0, lc};
Point(9) = {-delta, -delta, -0, lc};
Point(10) = {1 + delta, 1 + delta, -0, lc};
Point(11) = {-delta, 1 + delta, -0, lc};



Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line(6) = {11, 9};
Line(7) = {9, 8};
Line(8) = {8, 10};
Line(9) = {10, 11};

Spline(10) = {12, 14, 13, 15, 12};
Spline(11) = {17, 19, 18, 20, 17};

Line Loop(14) = {10};
Plane Surface(15) = {14};

Line Loop(20) = {11};
Plane Surface(20) = {20};

Line Loop(16) = {4, 1, 2, 3};
Plane Surface(17) = {16, 14, 20};
Line Loop(18) = {9, 6, 7, 8};
Plane Surface(19) = {16, 18};

//=============== LABELING ===============//
// Interface
Physical Line(12) = {10};
Physical Line(14) = {11};
//Physical Line(9) = {1, 2, 3, 4};
//Physical Line(13) = {6, 7, 8, 9};

// Omega_(...)
Physical Surface(1) = {15};
Physical Surface(2) = {20};
Physical Surface(3) = {17};
Physical Surface(4) = {19};
 

