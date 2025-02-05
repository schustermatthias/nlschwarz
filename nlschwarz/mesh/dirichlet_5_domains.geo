delta = 0.1;
cl__1 = 1;
lc = 1;
lc2 = 1;
lc3 = 1;
lc_fine = lc/2;
lc_point= lc/(15);


// Circle 1
Point(10) = {0.25, 0.05, 0, lc};
Point(11) = {0.25, 0.25, 0, lc};
Point(12) = {0.25, 0.45, 0, lc};

// Circle 2
Point(13) = {0.75, 0.05, 0, lc};
Point(14) = {0.75, 0.25, 0, lc};
Point(15) = {0.75, 0.45, 0, lc};

// Circle 3
Point(16) = {0.75, 0.55, 0, lc};
Point(17) = {0.75, 0.75, 0, lc};
Point(18) = {0.75, 0.95, 0, lc};

// Circle 4
Point(19) = {0.25, 0.55, 0, lc};
Point(20) = {0.25, 0.75, 0, lc};
Point(21) = {0.25, 0.95, 0, lc};

//------------------------------------------------------------------------------
// OMEGA_3
Point(1) = {0, 0, 0,  lc};
Point(2) = {0, 1, 0,  lc};
Point(3) = {1, 1, 0,  lc};
Point(4) = {1, 0, 0,  lc};

// Omega_I
Point(5) = {1. + delta, -delta, -0, lc};
Point(6) = {-delta, -delta, -0, lc};
Point(7) = {1 + delta, 1 + delta, -0, lc};
Point(8) = {-delta, 1 + delta, -0, lc};



Line(1) = {2, 1};
Line(2) = {1, 4};
Line(3) = {4, 3};
Line(4) = {3, 2};
Line(5) = {8, 6};
Line(6) = {6, 5};
Line(7) = {5, 7};
Line(8) = {7, 8};


Circle(10) = {10, 11, 12};
Circle(11) = {12, 11, 10};

Circle(12) =  {13, 14, 15};
Circle(13) =  {15, 14, 13};

Circle(14) =  {16, 17, 18};
Circle(15) =  {18, 17, 16};

Circle(16) =  {19, 20, 21};
Circle(17) =  {21, 20, 19};

Line Loop(20) = {10, 11};
Plane Surface(20) = {20};

Line Loop(21) = {12, 13};
Plane Surface(21) = {21};

Line Loop(22) = {14, 15};
Plane Surface(22) = {22};

Line Loop(23) = {16, 17};
Plane Surface(23) = {23};

Line Loop(24) = {4, 1, 2, 3};
Plane Surface(25) = {24, 20, 21, 22, 23};
Line Loop(26) = {8, 5, 6, 7};
Plane Surface(27) = {26, 24};

//=============== LABELING ===============//
// Interface
Physical Line(10) = {10, 11};
Physical Line(12) = {12, 13};
Physical Line(14) = {14, 15};
Physical Line(16) = {16, 17};
//Physical Line(9) = {1, 2, 3, 4};
//Physical Line(13) = {6, 7, 8, 9};

// Omega_(...)
Physical Surface(1) = {20};
Physical Surface(2) = {21};
Physical Surface(3) = {22};
Physical Surface(4) = {23};
Physical Surface(5) = {25};
Physical Surface(6) = {27};
 

