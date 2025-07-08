SetFactory("OpenCASCADE");

// –– mesh size parameters (coarse & uniform)
lc = 0.05;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// –– define two stacked boxes along z:
//    bottom half → retina (high absorption)
//    top half    → vitreous (low absorption)

// id 1: retina   (z ∈ [0.00, 0.50])
Box(1) = { 0, 0, 0.00,   1, 1, 0.50 };

// id 2: vitreous (z ∈ [0.50, 1.00])
Box(2) = { 0, 0, 0.50,   1, 1, 0.50 };

// –– tag each sub-volume for extraction in Feel++:
Physical Volume("retina")   = { 1 };
Physical Volume("vitreous") = { 2 };

// –– (Optional) tag the external cube faces if you need boundary markers:
Physical Surface("Left")   = { Surface{1,2}.FacesOfType("Xmin") };
Physical Surface("Right")  = { Surface{1,2}.FacesOfType("Xmax") };
Physical Surface("Front")  = { Surface{1,2}.FacesOfType("Ymin") };
Physical Surface("Back")   = { Surface{1,2}.FacesOfType("Ymax") };
Physical Surface("Bottom") = { Surface{1,2}.FacesOfType("Zmin") };
Physical Surface("Top")    = { Surface{1,2}.FacesOfType("Zmax") };
