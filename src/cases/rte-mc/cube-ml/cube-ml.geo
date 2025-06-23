SetFactory("OpenCASCADE");

// –– mesh size parameters (coarse & uniform)
lc = 0.05;
Mesh.CharacteristicLengthMin = lc;
Mesh.CharacteristicLengthMax = lc;

// –– define five stacked boxes along z (each thickness = 0.20):
//    retina  → z ∈ [0.00, 0.20]
//    vitreous→ z ∈ [0.20, 0.40]
//    lens    → z ∈ [0.40, 0.60]
//    aqueous → z ∈ [0.60, 0.80]
//    cornea  → z ∈ [0.80, 1.00]

// id 1: retina  (z ∈ [0.00, 0.20])
Box(1) = { 0, 0, 0.00,   1, 1, 0.20 };

// id 2: vitreous (z ∈ [0.20, 0.40])
Box(2) = { 0, 0, 0.20,   1, 1, 0.20 };

// id 3: lens     (z ∈ [0.40, 0.60])
Box(3) = { 0, 0, 0.40,   1, 1, 0.20 };

// id 4: aqueous  (z ∈ [0.60, 0.80])
Box(4) = { 0, 0, 0.60,   1, 1, 0.20 };

// id 5: cornea   (z ∈ [0.80, 1.00])
Box(5) = { 0, 0, 0.80,   1, 1, 0.20 };

// –– tag each sub-volume for extraction in Feel++:
Physical Volume("retina")   = { 1 };
Physical Volume("vitreous") = { 2 };
Physical Volume("lens")     = { 3 };
Physical Volume("aqueous")  = { 4 };
Physical Volume("cornea")   = { 5 };

// –– (Optional) tag the external cube faces if you need boundary markers:
Physical Surface("Left")   = { Surface{1,2,3,4,5}.FacesOfType("Xmin") };
Physical Surface("Right")  = { Surface{1,2,3,4,5}.FacesOfType("Xmax") };
Physical Surface("Front")  = { Surface{1,2,3,4,5}.FacesOfType("Ymin") };
Physical Surface("Back")   = { Surface{1,2,3,4,5}.FacesOfType("Ymax") };
Physical Surface("Bottom") = { Surface{1,2,3,4,5}.FacesOfType("Zmin") };
Physical Surface("Top")    = { Surface{1,2,3,4,5}.FacesOfType("Zmax") };
