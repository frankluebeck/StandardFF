#############################################################################
##
#A  StandardCyc.gd                                               Frank Lübeck
##  
##  The files StandardCyc.g{i,d} contain code to compute standard generators
##  of cyclic subgroups of multiplicative groups in standard finite fields.
##  


DeclareGlobalFunction("StdCycGen");

DeclareOperation("StandardCyclicGenerator", 
                             [IsStandardFiniteField, IsPosInt]);
DeclareAttribute("StandardPrimitiveRoot", IsStandardFiniteField);


