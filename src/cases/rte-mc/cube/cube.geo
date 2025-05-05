// cube.geo
SetFactory("OpenCASCADE");

// –– mesh size parameters (coarse & uniform)
lc = 0.5;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// –– build a 1×1×1 cube
Box(1) = { 0, 0, 0,   1, 1, 1 };

// –– tag the volume and its faces for boundary conditions or extraction
Physical Volume("Cube")    = { 1 };
Physical Surface("Left")   = { 1 };
Physical Surface("Right")  = { 2 };
Physical Surface("Front")  = { 3 };
Physical Surface("Back")   = { 4 };
Physical Surface("Bottom") = { 5 };
Physical Surface("Top")    = { 6 };