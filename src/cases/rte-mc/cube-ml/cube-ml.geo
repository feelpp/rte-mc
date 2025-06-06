// cube_layers.geo
SetFactory("OpenCASCADE");

// –– mesh size parameters (coarse & uniform)
lc = 0.1;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// –– define four stacked boxes along z:
//    bottom quarter → vitreous
//    next      → lens
//    next      → aqueous
//    top       → cornea

// id 1: vitreous (z ∈ [0.00, 0.25])
Box(1) = { 0, 0, 0.00,   1, 1, 0.25 };

// id 2: lens     (z ∈ [0.25, 0.50])
Box(2) = { 0, 0, 0.25,   1, 1, 0.25 };

// id 3: aqueous  (z ∈ [0.50, 0.75])
Box(3) = { 0, 0, 0.50,   1, 1, 0.25 };

// id 4: cornea   (z ∈ [0.75, 1.00])
Box(4) = { 0, 0, 0.75,   1, 1, 0.25 };

// –– tag each sub-volume for extraction in Feel++:
Physical Volume("vitreous") = { 1 };
Physical Volume("lens")     = { 2 };
Physical Volume("aqueous")  = { 3 };
Physical Volume("cornea")   = { 4 };

// –– (Optional) tag the external cube faces if you need boundary markers:
Physical Surface("Left")   = { Surface{1,2,3,4}.FacesOfType("Xmin")   };
Physical Surface("Right")  = { Surface{1,2,3,4}.FacesOfType("Xmax")   };
Physical Surface("Front")  = { Surface{1,2,3,4}.FacesOfType("Ymin")   };
Physical Surface("Back")   = { Surface{1,2,3,4}.FacesOfType("Ymax")   };
Physical Surface("Bottom") = { Surface{1,2,3,4}.FacesOfType("Zmin")   };
Physical Surface("Top")    = { Surface{1,2,3,4}.FacesOfType("Zmax")   };