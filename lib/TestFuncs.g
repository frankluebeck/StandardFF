#############################################################################
##
#A  TestFuncs.g                                                  Frank Lübeck
##  
##  The file TestFuncs.g contains various functions for testing the package.
##  


##  <#GAPDoc Label="TestLoops">
##  <ManSection>
##  <Heading>Computing all fields in various ranges</Heading>
##  <Func Name="AllPrimeDegreePolynomials" Arg="p, bound"/>
##  <Func Name="AllFF" Arg="p, bound"/>
##  <Func Name="AllPrimitiveRoots" Arg="p, bound"/>
##  <Func Name="AllPrimitiveRootsCANFACT" Arg=""/>
##  <Func Name="AllFieldsWithConwayPolynomial" Arg='["ConwayGen",]["MiPo"]'/>
##  <Description>
##  These function  compute all fields  in some range, sometimes  with further
##  data. All functions  return a list with some timings  and print a log-file
##  in the current directory.
##  <P/> 
##  <Ref    Func="AllPrimeDegreePolynomials"/>   computes    all   irreducible
##  polynomials  of prime  degree needed  for the  construction of  all finite
##  fields of order <M><A>p</A>^i</M>, <M>1 \leq i \leq <A>bound</A></M>. This
##  is the most time consuming part in the construction of the fields.
##  <P/> 
##  <Ref  Func="AllFF"/> computes  all  <C>FF(p,i)</C> for  <M>1  \leq i  \leq
##  <A>bound</A></M>. When  the previous  function was  called before  for the
##  same range, this function spends most of its time by computing the minimal
##  polynomials of the standardized primitive elements of <C>FF(p,i)</C>.
##  <P/>
##  <Ref Func="AllPrimitiveRoots"/> computes  the standardized primitive roots
##  in <C>FF(p,i)</C>  for <M>1  \leq i \leq  <A>bound</A></M>. The  most time
##  consuming cases  are when a  large prime divisor <M>r</M>  of <M>p^i-1</M>
##  already divides <M>p^j-1</M>  for some <M>j &lt; i</M>  (but then <M>r</M>
##  divides <M>i/j</M>). Cases where &GAP; cannot factorize <M>p^i-1</M> (that
##  is <M>i</M> is not contained in <C>CANFACT[p]</C>) are skipped.
##  <P/>
##  <Ref  Func="AllPrimitiveRootsCANFACT"/>  does  the same  as  the  previous
##  function for all pairs <M>p, i</M> stored in <Ref Var="CANFACT"/>.
##  <P/>
##  <Ref  Func="AllFieldsWithConwayPolynomial"/>  computes all  <C>FF(p,i)</C>
##  for     the     cases     where    &GAP;     knows     the     precomputed
##  <C>ConwayPolynomial(p,i)</C>.      With     the      optional     argument
##  <C>"ConwayGen</C>  the   function  computes   for  all  fields   the  <Ref
##  Attr="SteinitzPairConwayGenerator"/>   and   writes   it   into   a   file
##  <F>SteinitzPairConway</F>.  With   the  optional   argument  <C>"MiPo"</C>
##  the   function   also   computes    the   minimal   polynomials   of   the
##  <Ref   Attr="StandardPrimitiveRoot"/>   and   writes   it   to    a   file
##  <F>MiPoPrimitiveRoots</F> (these  polynomials have the  same compatibility
##  properties as Conway polynomials).
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
AllPrimeDegreePolynomials := function(p, bound, cache...)
  local nam, out, T, times, rr, k, t, pl, tt, r, outc, snfpd;
  if "NTL" in cache then
    snfpd := SteinitzNumberForPrimeDegreeNTL;
  else
    snfpd := SteinitzNumberForPrimeDegree;
  fi;
  if Length(cache) > 0 and cache[1] = true then
    cache := true;
  else
    cache := false;
  fi;
  nam := Concatenation("p",String(p),"rk_",String(bound),".log");
  out := OutputTextFile(nam, false);
  if cache then
    nam := Concatenation("p_",String(p),"_stpp.g");
    outc := OutputTextFile(nam, true);
    PrintTo(outc, "\nif not IsBound(STPP) then STPP := []; fi;\n");
  fi;
  
  T := Runtime();
  times := [];
  for r in [1..bound] do
    if IsPrime(r) then
      rr := r;
      k := 1;
      while rr <= bound do
        Print([p,r,k],"\c");
        t := Runtime();
        PrintTo(out,"# ",p,", ",r,", ",k,": \c");
        #pl := StandardPrimeDegreePolynomial(p,r,k);
        pl := snfpd(p,r,k);
        if cache then
          PrintTo(outc, "Add(STPP,",[p,r,k,pl],");\n");
        fi;
        tt := Runtime()-t;
        PrintTo(out, tt, "\n");
        Add(times, [p,r,k,tt]);
        rr := rr*r;
        k := k+1;
      od;
    fi;
  od;
  PrintTo(out, "\n\n# Total time: ",StringTime(Runtime()-T), "\n\n");
  PrintTo(out, "times := \n", times, ";\n\n");
  CloseStream(out);
  if cache then
    CloseStream(outc);
  fi;
  return times;
end;
AllFF := function(p, bound)
  local nam, out, T, times, t, pl, tt, n;
  nam := Concatenation("p",String(p),"nFF_",String(bound),".log");
  out := OutputTextFile(nam, false);
  T := Runtime();
  times := [];
  for n in [1..bound] do
    Print([p,n],"\c");
    t := Runtime();
    PrintTo(out, "# FF ",p,", ",n,": \c");
    pl := FF(p,n);
    tt := Runtime()-t;
    PrintTo(out, tt, "\n");
    Add(times, [p,n,tt]);
  od;
  PrintTo(out, "\n\n# Total time: ",StringTime(Runtime()-T), "\n\n");
  PrintTo(out, "times := \n", times, ";\n\n");
  CloseStream(out);
  return times;
end;
AllPrimitiveRoots := function(p, bound)
  local nam, out, T, times, t, pl, tt, n;
  nam := Concatenation("p",String(p),"nPrimRoot_",String(bound),".log");
  out := OutputTextFile(nam, false);
  T := Runtime();
  times := [];
  for n in [1..bound] do
    if n in CANFACT[p] then
      Print([p,n],"\c");
      t := Runtime();
      PrintTo(out, "# FF ",p,", ",n,": \c");
      pl := StandardPrimitiveRoot(FF(p,n));
      tt := Runtime()-t;
      PrintTo(out, tt, "\n");
      Add(times, [p,n,tt]);
    fi;
  od;
  PrintTo(out, "\n\n# Total time: ",StringTime(Runtime()-T), "\n\n");
  PrintTo(out, "times := \n", times, ";\n\n");
  CloseStream(out);
  return times;
end;

AllPrimitiveRootsCANFACT := function()
  local pr, nam, out, T, times, t, pl, tt, p, n;
  pr := Filtered([1..Length(CANFACT)], i-> IsBound(CANFACT[i]));
  nam := Concatenation("PrimRootCANFACT.log","");
  out := OutputTextFile(nam, false);
  T := Runtime();
  times := [];
  for p in pr do
    for n in CANFACT[p]  do
      Print([p,n],"\c");
      t := Runtime();
      PrintTo(out, "# FF ",p,", ",n,": \c");
      pl := StandardPrimitiveRoot(FF(p,n));
      tt := Runtime()-t;
      PrintTo(out, tt, "\n");
      Add(times, [p,n,tt]);
    od;
  od;
  PrintTo(out, "\n\n# Total time: ",StringTime(Runtime()-T), "\n\n");
  PrintTo(out, "times := \n", times, ";\n\n");
  CloseStream(out);
  return times;
end;
AllFieldsWithConwayPolynomial := function(arg...)
  local nam, out, T, times, cd, t, f, tt, p, j, x, st, pol;
  if "MiPo" in arg then
    PrintTo("MiPoPrimitiveRoots", "mp := [];\n");
  fi;
  if "ConwayGen" in arg then
    PrintTo("SteinitzPairsConway", "li := [];\n");
  fi;
  nam := Concatenation("PrimRootConway.log","");
  out := OutputTextFile(nam, false);
  T := Runtime();
  times := [];
  for p in [2..10000] do
    if IsPrimeInt(p) then
      ConwayPol(p,2);
      cd := CONWAYPOLDATA[p];
      for j in [1..Length(cd)] do
        if IsBound(cd[j]) then
          Print([p,j],"\c");
          t := Runtime();
          PrintTo(out, "# FF ",p,", ",j,": \c");
          f := FF(p, j);
          x := StandardPrimitiveRoot(f); 
          if "MiPo" in arg then
            pol := MinimalPolynomialByBerlekampMasseyShoup(x);
            st :=
              ValuePol(List(CoefficientsOfUnivariatePolynomial(pol),IntFFE),p);
            AppendTo("MiPoPrimitiveRoots", "Add(mp, ", [p,j,st],");\n");
          fi;
          if "ConwayGen" in arg then
            st := SteinitzPairConwayGenerator(f);
            AppendTo("SteinitzPairsConway", "Add(li, ", [p,st],");\n");
          fi;
          tt := Runtime()-t;
          PrintTo(out, tt, "\n");
          Add(times, [p,j,tt]);
        fi;
      od;
    fi;
  od;
  PrintTo(out, "\n\n# Total time");
  if "ConwayGen" in arg then
    PrintTo(out," (with Conway generator)");
  fi;
  PrintTo(out, ": ",StringTime(Runtime()-T), "\n\n");
  PrintTo(out, "times := \n", times, ";\n\n");
  CloseStream(out);
  return times;
end;


