#############################################################################
##
#A  read.g              StandardFF package                       Frank Lübeck
##
##  Reading the library of the package.
##  

ReadPackage( "StandardFF", "lib/TowerArith.gi");
ReadPackage( "StandardFF", "lib/Util.gi");
ReadPackage( "StandardFF", "lib/IsIrred.gi");
ReadPackage( "StandardFF", "lib/StandardFF.gi");
ReadPackage( "StandardFF", "lib/StandardCyc.gi");
ReadPackage( "StandardFF", "lib/CANFACT.g");

ReadPackage( "StandardFF", "lib/TestFuncs.g");

# redefine the Info handler and output, we use the plain one from GAPDoc
SetInfoHandler(InfoStandardFF, PlainInfoHandler);
SetInfoOutput(InfoStandardFF, "*errout*");

# user preference
DeclareUserPreference(rec(
  name := "StandardFFUseCache",
  description := ["""
If set to 'true' the package StandardFF will use cached data from its \\
'data' directory. Otherwise it will compute everything from scratch. \\
The default is 'false'.
"""],
  default := false,
  values := [true, false],
  multi := false,
  package := "StandardFF"));

InstallAtExit( function()
  local fn, c;
  if UserPreference("StandardFFUseCache") = true then
    for c in ["IRRPRIMDEGCACHE", "PRIMROOTCACHE", "MIPOPRIMCACHE"] do
      if IsBound(SFFHelper.(c)) then
        fn := Filename(DirectoriesPackageLibrary("StandardFF","data"), c);
        if fn <> fail then
          PrintTo(fn, "SFFHelper.",c," := \n",
                      SFFHelper.(c),
                      ";\n\n");
        fi;
      fi;
    od;
  fi;
end);
