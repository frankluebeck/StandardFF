#############################################################################
##
#A  TowerArith.gd                                                Frank Lübeck
##
##  The  files TowerArith.g{i,d}  contain code  to generate  and compute  in
##  towers of finite fields.
##  

# for multivariate polynomial wrt one indeterminate
DeclareOperation("Degree", [IsPolynomial, IsObject]);
DeclareOperation("Indets", [IsObject]);

DeclareCategory( "IsFiniteFieldTower", 
                  IsField and IsFinite and IsAttributeStoringRep);

DeclareCategory( "IsFiniteFieldTowerElement", 
                          IsRingElementWithInverse and 
                          IsAssociativeElement and
                          IsAdditivelyCommutativeElement and
                          IsCommutativeElement and
                          IsAttributeStoringRep);

BindGlobal("FiniteFieldTowerElementFamily", 
     NewFamily("FiniteFieldTowerElementFamily", IsFiniteFieldTowerElement));
Setter(IsUFDFamily)(FiniteFieldTowerElementFamily, true);
BindGlobal("FiniteFieldTowerFamily", 
     CollectionsFamily(FiniteFieldTowerElementFamily));
BindGlobal("FiniteFieldTowerType", 
     NewType(FiniteFieldTowerFamily, IsFiniteFieldTower));
BindGlobal("FiniteFieldTowerElementType", 
     NewType(FiniteFieldTowerElementFamily, IsFiniteFieldTowerElement));

DeclareGlobalFunction("FiniteFieldTower");
DeclareGlobalFunction("ExtendFiniteFieldTower");
DeclareGlobalFunction("PreviousFiniteFieldTower");
DeclareAttribute("BaseFiniteFieldTower", IsFiniteFieldTower);

DeclareOperation("FiniteFieldTowerElement", [IsFiniteFieldTower, IsPolynomial]);
DeclareOperation("FiniteFieldTowerElement", [IsFiniteFieldTower, IsList]);
DeclareOperation("FiniteFieldTowerElement", 
                               [IsFiniteFieldTower, IsFiniteFieldTowerElement]);
DeclareGlobalFunction("FiniteFieldTowerElementVector");
DeclareGlobalFunction("PolynomialFiniteFieldTowerElement");

DeclareOperation("Reduce", [IsFiniteFieldTowerElement]);

DeclareAttribute("ReducedGeneratorPowers", IsFiniteFieldTower);

DeclareOperation("AsVector", [IsFiniteFieldTowerElement]);

# use only for towers with prime field as base
DeclareOperation("SteinitzNumber", [IsFiniteFieldTowerElement]);
DeclareOperation("SteinitzNumber", [IsFFE]);
DeclareOperation("ElementSteinitzNumber", [IsFiniteFieldTower, IsInt]);

# basis of monomials corresponding to tower
DeclareAttribute("TowerBasis", IsFiniteFieldTower);
