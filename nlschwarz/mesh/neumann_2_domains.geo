delta = 0.1;
cl__1 = 1;
lc = 1;

//----------------------------------------
// Omega
Point(1) = {0, 0, 0, 1};
Point(2) = {1, 0, 0, 1};
Point(3) = {1, 1, 0, 1};
Point(4) = {0, 1, 0, 1};

// Omega_I
Point(8) = {-delta, 0, 0, lc};
Point(9) = {0, -delta, 0, lc};
Point(10) = {1, -delta, 0, lc};
Point(11) = {1 + delta, 0, 0, lc};

Point(14) = {1 + delta, 1, 0, lc};
Point(15) = {1, 1 + delta, 0, lc};
Point(16) = {0, 1 + delta, 0, lc};
Point(17) = {-delta, 1, 0, lc};

// for midline
Point(5)  = {0.5,       0.0, 0, 1};
Point(6)  = {0.5,       1.0, 0, 1};
Point(12) = {0.5,    -delta, 0, 1};
Point(13) = {0.5, 1 + delta, 0, 1};

//Lines around Omega(_1 and _2)
Line(1) = {1, 5};
Line(2) = {5, 2};
Line(3) = {2, 3};
Line(4) = {3, 6};
Line(5) = {6, 4};
Line(6) = {4, 1};

//Lines around Omega_I
Circle(21) = {8, 1, 9};
Line(22) = {9, 12};
Line(23) = {12, 10};
Circle(24) = {10, 2, 11};
Line(25) = {11, 14};
Circle(26) = {14, 3, 15};
Line(27) = {15, 13};
Line(28) = {13, 16};
Circle(29) = {16, 4, 17};
Line(30) = {17, 8};

//midlines
Line(13) = {5, 6};
Line(14) = {12, 5};
Line(15) = {6, 13};

//Omega_1
Curve Loop(1) = {1, 13, 5, 6};
Plane Surface(1) = {1};

//Omega_2
Curve Loop(3) = {2, 3, 4, -13};
Plane Surface(3) = {3};

//Gamma_1(with tag 2)
Curve Loop(2) = {22, 14, -1, -6, -5, 15, 28, 29, 30, 21};
Plane Surface(2) = {2};

//Gamma_2(with tag 4)
Curve Loop(4) = {23, 24, 25, 26, 27, -15, -4, -3, -2, -14};
Plane Surface(4) = {4};

//Tag every surface
Physical Surface(1) = {1};
Physical Surface(2) = {2};
Physical Surface(3) = {3};
Physical Surface(4) = {4};

Physical Curve(1) = {13};
Physical Curve(2) = {14};
Physical Curve(3) = {15};

