##############################################################################
##
#A  generateCANFACT.g                                             Frank Lübeck
##
##  This file can be read by GAP to produce an update of the list 'CANFACT' in
##  'CANFACT.g'. For example after an update of the list of Brent factors, see
##  '?FetchMoreFactors'.
##  
##  The idea is to fork GAP instances which try for 5 seconds to factorize p^n
##  - 1. If successful the pair (p,n) is included in CANFACT.
##


SetInfoLevel( InfoPrimeInt, 0 );
out := OutputTextFile("CANFACT", false);
SetPrintFormattingStatus(out, false);
AppendTo(out, "CANFACT := [];\n");
for p in Filtered([1..10000], IsPrime)  do
  Print(p," [\c");
  if p < 10 then
    max := 2000;
  elif p < 100 then 
    max := 500;
  else 
    max := 100;
  fi;
  AppendTo(out, "CANFACT[",p,"]:=[];\n");
  l := [];
  for n in [1..max] do
    Print([p,n]," \c");
    pf := ParTakeFirstResultByFork([Factors],[[p^n-1]],
          rec(TimeOut:=rec(tv_sec:=5,tv_usec:=0)));
    if Length(pf) > 0 then
      AppendTo(out, "Add(CANFACT[",p,"],",n,");\n");
    fi;
  od;
od;
CloseStream(out);
Print("""Read the file 'CANFACT' and print the variable CANFACT to the file
    CANFACT.g
to store the data more compactly.
""");



