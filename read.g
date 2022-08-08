#############################################################################
##
#A  read.g              StandardFF package                       Frank Lübeck
##
##  Reading the library of the package.
##  

# workaround for some (not very realistic) CI tests
if not IsBound(BRENTFACTORS) then
  BRENTFACTORS := 0;
  WriteBrentFactorsFiles := 0;
fi;
if not IsBound(IO_PipeThrough) then
  IO_PipeThrough := 0;
fi;


ReadPackage( "StandardFF", "lib/Util.gi");
ReadPackage( "StandardFF", "lib/IsIrred.gi");
ReadPackage( "StandardFF", "lib/StandardFF.gi");
ReadPackage( "StandardFF", "lib/StandardCyc.gi");
ReadPackage( "StandardFF", "lib/CANFACT.g");

ReadPackage( "StandardFF", "lib/TestFuncs.g");

# redefine the Info handler and output, we use the plain one from GAPDoc
SetInfoHandler(InfoStandardFF, PlainInfoHandler);
SetInfoOutput(InfoStandardFF, "*errout*");


# part of workaround above:
if BRENTFACTORS = 0 then
  Unbind(BRENTFACTORS);
  Unbind(WriteBrentFactorsFiles);
fi;
if IO_PipeThrough = 0 then
  Unbind(IO_PipeThrough);
fi;


