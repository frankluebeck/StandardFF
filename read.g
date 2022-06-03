#############################################################################
##
#A  read.g              StandardFF package                       Frank Lübeck
##
##  Reading the library of the package.
##  

ReadPackage( "StandardFF", "lib/Util.gi");
ReadPackage( "StandardFF", "lib/IsIrred.gi");
ReadPackage( "StandardFF", "lib/StandardFF.gi");
ReadPackage( "StandardFF", "lib/StandardCyc.gi");
ReadPackage( "StandardFF", "lib/CANFACT.g");

ReadPackage( "StandardFF", "lib/TestFuncs.g");

# redefine the Info handler and output, we use the plain one from GAPDoc
SetInfoHandler(InfoStandardFF, PlainInfoHandler);
SetInfoOutput(InfoStandardFF, "*errout*");

