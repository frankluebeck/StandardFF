#############################################################################
##
#A  StandardFF.gd                                                Frank Lübeck
##  
##  The files StandardFF.g{i,d}  contain code to compute  standard models of
##  finite fields.  Fields are  returned as field  towers with  prime degree
##  indices, and as single extensions over the prime field.
##  

# standard finite fields as towers of prime degree extensions
DeclareGlobalFunction("StandardPrimeDegreePolynomial");
DeclareGlobalFunction("StandardFiniteFieldTower");
# standard finite fields as simple extensions over prime field GF(p)
DeclareGlobalFunction("StandardFiniteFieldNonSparse");
DeclareGlobalFunction("StandardFiniteField");
DeclareSynonym("FF", StandardFiniteField);

DeclareCategory("IsStandardFiniteFieldTower", IsFiniteFieldTower);
DeclareProperty("IsStandardFiniteField", 
                         IsField and IsFinite and IsAlgebraicExtension);
DeclareCategory("IsStandardFiniteFieldElement", IsAlgebraicElement);

DeclareAttribute("Tower", IsStandardFiniteField);
DeclareAttribute("PrimitivePowersInTowerBasis", IsStandardFiniteField);
DeclareAttribute("TowerBasis", IsStandardFiniteField);
DeclareAttribute("TowerBasisMap", IsStandardFiniteField);


DeclareOperation("ToTowerElement", [IsStandardFiniteField, IsObject]);
DeclareOperation("FromTowerElement", [IsStandardFiniteField, IsObject]);
DeclareOperation("SteinitzNumber", [IsStandardFiniteField, IsRingElement]);
DeclareOperation("ElementSteinitzNumber", [IsStandardFiniteField, IsInt]);

DeclareOperation("AsVector", [IsStandardFiniteFieldElement]);
DeclareAttribute("AsPolynomial", IsStandardFiniteFieldElement);
DeclareOperation("ElementVector", [IsStandardFiniteField, IsRowVector]);
DeclareOperation("ElementPolynomial", [IsStandardFiniteField, IsPolynomial]);

# records for helper functions
SFFHelper := rec();
MakeReadOnlyGlobal("SFFHelper");
TowerMon := rec();
MakeReadOnlyGlobal("TowerMon");

# embeddings
DeclareGlobalFunction("StdMon");
DeclareGlobalFunction("StandardIsomorphismGF");
DeclareCategory("IsStandardFiniteFieldEmbedding", IsFieldHomomorphism);
DeclareOperation("SteinitzPair", [IsStandardFiniteFieldElement]);
# use only for standard finite field towers
DeclareOperation("SteinitzPair", [IsFiniteFieldTowerElement]);


DeclareAttribute("SteinitzPairConwayGenerator", IsStandardFiniteField);

# helper to identify GF(p)'s in method selection
DeclareProperty("IsStandardPrimeField", IsField and IsFinite);
