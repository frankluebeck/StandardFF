#############################################################################
##
#A  StandardFF.gd                                                Frank Lübeck
##  
##  The files StandardFF.g{i,d}  contain code to compute  standard models of
##  finite fields.  Fields are  returned as field  towers with  prime degree
##  indices, and as single extensions over the prime field.
##  

# standard finite fields as towers of prime degree extensions
DeclareGlobalFunction("SteinitzNumberForPrimeDegree");
DeclareGlobalFunction("SteinitzNumberForPrimeDegreeNTL");
DeclareGlobalFunction("StandardPrimeDegreePolynomial");
DeclareGlobalFunction("_ExtensionWithTowerBasis");
# standard finite fields as simple extensions over prime field GF(p)
DeclareGlobalFunction("StandardFiniteField");
DeclareSynonym("FF", StandardFiniteField);

DeclareProperty("IsStandardFiniteField", 
                         IsField and IsFinite and IsAlgebraicExtension);
DeclareCategory("IsStandardFiniteFieldElement", 
                         IsAlgebraicElement and IsKroneckerConstRep and
                         IsNoImmediateMethodsObject);

DeclareAttribute("PrimitivePowersInTowerBasis", IsStandardFiniteField);
DeclareAttribute("TowerBasis", IsStandardFiniteField);

DeclareOperation("SteinitzNumber", [IsStandardFiniteField, IsRingElement]);
DeclareOperation("ElementSteinitzNumber", [IsStandardFiniteField, IsInt]);

DeclareAttribute("AsVector", IsStandardFiniteFieldElement);
DeclareAttribute("GeneratorMonomials", IsStandardFiniteField);
DeclareAttribute("TowerBasisMonomials", IsStandardFiniteField);
DeclareAttribute("AsPolynomial", IsStandardFiniteFieldElement);
DeclareOperation("ElementVector", [IsStandardFiniteField, IsRowVector]);
DeclareOperation("ElementPolynomial", [IsStandardFiniteField, IsPolynomial]);

# embeddings
DeclareGlobalFunction("StdMon");
DeclareGlobalFunction("StdMonDegs");
DeclareGlobalFunction("StdMonMap");
DeclareGlobalFunction("EmbedSteinitz");
DeclareCategory("IsStandardFiniteFieldEmbedding", IsFieldHomomorphism);
DeclareOperation("SteinitzPair", [IsStandardFiniteFieldElement]);

DeclareAttribute("SteinitzPairConwayGenerator", IsStandardFiniteField);
DeclareGlobalFunction("StandardIsomorphismGF");

# helper to identify GF(p)'s in method selection
DeclareProperty("IsStandardPrimeField", IsField and IsFinite);

# creating elements
DeclareOperation("ZZ", [IsInt, IsInt, IsList]);

# arithmetic in algebraic closure
DeclareGlobalFunction("MoveToCommonStandardField");
DeclareGlobalFunction("MoveToSmallestStandardField");
