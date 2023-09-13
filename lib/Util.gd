#############################################################################
##
#A  SimpleRandom.gd                                              Frank Lübeck
##  
##  The files Util.g{i,d}  contain code for some utiltiy functions.
##  

DeclareGlobalFunction("StandardAffineShift");

DeclareGlobalFunction("FindLinearCombination");

DeclareGlobalFunction("OrderModBound");

DeclareOperation("InvModCoeffs", [IsList, IsList]);

if not IsBound(DLog) then
  DeclareGlobalName("DLogShanks");
  DeclareGlobalName("DLog");
fi;

DeclareGlobalFunction("FindConjugateZeroes");
DeclareGlobalFunction("FindConjugateZeroesChar2");
DeclareGlobalFunction("ZeroesConway");

DeclareGlobalFunction("BerlekampMassey");
DeclareGlobalFunction("MinimalPolynomialByBerlekampMassey");
DeclareGlobalFunction("MinimalPolynomialByBerlekampMasseyShoup");

DeclareGlobalFunction("FetchMoreFactors");

DeclareGlobalFunction("StandardValuesBrauerCharacter");
DeclareGlobalFunction("IsGaloisInvariant");
DeclareGlobalFunction("SmallestDegreeFrobeniusCharacterValue");
DeclareGlobalFunction("StandardFrobeniusCharacterValue");

# for multivariate polynomial wrt one indeterminate
DeclareOperation("Degree", [IsPolynomial, IsObject]);
DeclareOperation("Indets", [IsObject]);

