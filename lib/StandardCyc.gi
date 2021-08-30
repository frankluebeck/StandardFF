#############################################################################
##
#A  StandardCyc.gi                                               Frank Lübeck
##  
##  The files StandardCyc.g{i,d} contain code to compute standard generators
##  of cyclic subgroups of multiplicative groups in standard finite fields.
##  


# p, r: prime numbers
# bound: multiple of order of p mod r^k
# returns standard generator of order r^k in \bar GF(p) as Steinitz pair
InstallGlobalFunction(StdCycGen, function(p, r, k, bound)
  local F, cgens, res, l, m, mm, t, K, count, st, a, am, pl, b, 
        z, bb, aa, e, c, j, min, nn, s;
  # we cache the results
  F := GF(p);
  if not IsBound(F!.prcycgens) then
    F!.prcycgens := [];
  fi;
  cgens := F!.prcycgens;
  res := First(cgens, a-> a[1]=r and a[2]=k);
  if res <> fail then
    return res[3];
  fi;
  # smallest degree containing r-elements
  l := OrderModBound(p, r, bound);
  m := (p^l-1)/r;
  mm := m;
  t := 1;
  while m mod r = 0 do
    m := m/r;
    t := t+1;
  od;
  if t = k then
    # easy case:
    # we search for element whose m-th power has order p^t,
    # testing pseudo random elements  
    K := StandardFiniteField(p, l);
    # we use standard affine shift for pseudo random sequence
    count := 1;
    st := StandardAffineShift(Size(K), count);
    a := ElementSteinitzNumber(K, st);
    am := a^m;
    while IsZero(am) or IsOne(am^(r^(t-1))) do
      Info(InfoStandardFF, 2, "!\c");      
      count := count+1;
      st := StandardAffineShift(Size(K), count);
      a := ElementSteinitzNumber(K, st);
      am := a^m;
    od;
    pl := SteinitzPair(K, am);
    Add(cgens, [r, k, pl]);
    return pl;
  elif t > k then
    # just embed from smaller field
    K := StandardFiniteField(p, l);
    a := ElementSteinitzNumber(K, 
              SteinitzNumber(K, StdCycGen(p, r, t, bound)));
    am := a^(r^(t-k));
    pl := SteinitzPair(K, am);
    Add(cgens, [r, k, pl]);
    return pl;
  else
    # expensive case:
    # for r odd or r=2 and p = 1 (4):
    # have to compute an r-th root of element of order r^(k-1)
    # exception if r=2 and p = 3 (4):take recursively smaller square root,
    # starting with -1 in GF(p).
    l := OrderModBound(p, r^k, bound);
    K := StandardFiniteField(p, l);
    # this element generates a subfield of index r
    b := ElementSteinitzNumber(K, 
                    SteinitzNumber(K, StdCycGen(p, r, k-1, bound)));
    # find random element of order r^k
    m := (Size(K)-1)/r^k;
    a := Random(K)^m;
    while IsZero(a) or IsOne(a^(r^(k-1))) do
      a := Random(K)^m;
    od;
    # r-th root of one
    z := a^(r^(k-1));
    # find e with a^e is r-th root of b
    # (the j below are the coefficients of e in base r)
    bb := b;
    aa := a^r;
    e := 0;
    for s in [0..k-2] do
      c := bb^(r^(k-s-2));
      j := DLogShanks(z, c, r);
      e := e + j*r^s;
      bb := bb / aa^(j*r^s);
    od;
    a := a^e;
    # now a is an r-th root of b, we choose the one with smallest Steinitz
    # number
    aa := a;
    min := [aa, SteinitzNumber(K, aa)];
    for j in [1..r-1] do
      aa := aa*z;
      nn := SteinitzNumber(K, aa);
      if nn < min[2] then
        min := [aa, nn];
      fi;
    od;
    pl := SteinitzPair(K, min[1]);
    Add(cgens, [r, k, pl]);
    return pl;
  fi;
end);

##  <#GAPDoc Label="StandardCyclicGenerator">
##  <ManSection>
##  <Oper Name="StandardCyclicGenerator" Arg="F[, m]"/>
##  <Attr Name="StandardPrimitiveRoot" Arg="F"/>
##  <Returns>an element of <A>F</A> or <K>fail</K> </Returns>
##  <Description>
##  The  argument  <A>F</A>  must  be   a  standard  finite  field  (see  <Ref
##  Func="FF"/>) and <A>m</A> a positive  integer. If <A>m</A> does not divide
##  <M>|<A>F</A>|  -  1</M>  the  function returns  <K>fail</K>.  Otherwise  a
##  standardized element <M>x_{<A>m</A>}</M> of order <A>m</A> is returned, as
##  described above.
##  <P/> 
##  The  argument <A>m</A>  is optional,  if not  given its  default value  is
##  <M>|<A>F</A>|  -  1</M>. In  this  case  <M>x_{<A>m</A>}</M> can  also  be
##  computed with the attribute <Ref Attr="StandardPrimitiveRoot"/>.
##  <Example>gap> F := FF(67, 18); # Conway polynomial was hard to compute
##  FF(67, 18)
##  gap> x := PrimitiveElement(F);
##  x18
##  gap> xprim := StandardPrimitiveRoot(F);;
##  gap> k := (Size(F)-1) / Order(x);
##  6853662165340556076084083497526
##  gap> xm := StandardCyclicGenerator(F, Order(x));;
##  gap> xm = xprim^k;
##  true
##  gap> d := DLog(xm, x, Order(x));
##  31
##  gap> xm^d = x;
##  true
##  gap> F := FF(23, 201); # factorization of (|F| - 1) not known 
##  FF(23, 201)
##  gap> m:=79*269*67939;
##  1443771689
##  gap> (Size(F)-1) mod m;
##  0
##  gap> OrderMod(23, m);
##  201
##  gap> xm := StandardCyclicGenerator(F, m);;
##  gap> IsOne(xm^m);
##  true
##  gap> ForAll(Factors(m), r-> not IsOne(xm^(m/r)));
##  true
##  gap> F := FF(7,48);
##  FF(7, 48)
##  gap> K := FF(7,12);
##  FF(7, 12)
##  gap> emb := Embedding(K, F);;
##  gap> x := StandardPrimitiveRoot(F);;
##  gap> y := StandardPrimitiveRoot(K);;
##  gap> y^emb = x^((Size(F)-1)/(Size(K)-1));
##  true
##  gap> v := Indeterminate(FF(7,1), "v");
##  v
##  gap> px := MinimalPolynomial(FF(7,1), x, v);;
##  gap> py := MinimalPolynomial(FF(7,1), y, v);
##  v^12+Z(7)*v^11+Z(7)^2*v^10-v^9+v^6-v^4-v^3+Z(7)^5*v^2+Z(7)*v+Z(7)^5
##  gap> Value(py, PowerMod(v, (Size(F)-1)/(Size(K)-1), px)) mod px;
##  0*Z(7)
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
GAPInfo.tmpfun := function(K, ord)
  local p, n, fac, ppgens, res, num;
  if ord = 1 then
    return One(K);
  fi;
  p := Characteristic(K);
  n := DegreeOverPrimeField(K);
  if (p^n-1) mod ord <> 0 then
    return fail;
  fi;
  fac := Collected(Factors(ord));
  # standard generators of prime power order
  ppgens := List(fac, a-> StdCycGen(p, a[1], a[2], n));
  ppgens := List(ppgens, a-> SteinitzNumber(K, a));
  ppgens := List(ppgens, a-> ElementSteinitzNumber(K, a));

  # product has correct order, we take power to find element corresponding
  # to 1/ord
  res := Product(ppgens);
  num := Sum(fac, a-> 1/a[1]^a[2]) * ord;
  res := res^(1/num mod ord);
  return res;
end;
InstallMethod(StandardCyclicGenerator, 
     ["IsStandardFiniteField", "IsPosInt"], GAPInfo.tmpfun);
InstallOtherMethod(StandardCyclicGenerator, 
     ["IsStandardPrimeField", "IsPosInt"], GAPInfo.tmpfun);


GAPInfo.tmpfun := function(K)
  local p, n, res;
  if UserPreference("StandardFFUseCache") <> true then
    return StandardCyclicGenerator(K, Size(K)-1);
  fi;
  if not IsBound(SFFHelper.PRIMROOTCACHE) then
    ReadPackage("StandardFF", "data/PRIMROOTCACHE");
  fi;
  p := Characteristic(K);
  n := DegreeOverPrimeField(K);
  if IsBound(SFFHelper.PRIMROOTCACHE[p]) and 
                 IsBound(SFFHelper.PRIMROOTCACHE[p][n]) then
    return ElementSteinitzNumber(K, SFFHelper.PRIMROOTCACHE[p][n]);
  fi;
  res := StandardCyclicGenerator(K, Size(K)-1);
  if not IsBound(SFFHelper.PRIMROOTCACHE[p]) then
    SFFHelper.PRIMROOTCACHE[p] := [];
  fi;
  SFFHelper.PRIMROOTCACHE[p][n] := SteinitzNumber(res);
  return res;
end;
InstallMethod(StandardPrimitiveRoot, ["IsStandardFiniteField"],
GAPInfo.tmpfun);
InstallOtherMethod(StandardPrimitiveRoot, ["IsStandardPrimeField"],
GAPInfo.tmpfun);
Unbind(GAPInfo.tmpfun);
InstallMethod(PrimitiveRoot, ["IsStandardFiniteField"], StandardPrimitiveRoot);


InstallOtherMethod(StandardCyclicGenerator, ["IsStandardFiniteField"],
  StandardPrimitiveRoot);

# Log with respect to standard primitive root
InstallOtherMethod(Log, ["IsStandardFiniteFieldElement"], function(x)
  local F;
  F := FamilyObj(x)!.wholeField;
  return DLog(StandardPrimitiveRoot(F), x, Factors(Size(F)-1));
end);
# using order of F^\times
InstallOtherMethod(Order, ["IsStandardFiniteFieldElement"], function(x)
  local F, m, fac, res, k, xx, a;
  F := FamilyObj(x)!.wholeField;
  m := Size(F)-1;
  fac := Collected(Factors(m));
  res := m;
  for a in fac do
    k := a[2];
    xx := x^(m/a[1]^a[2]);
    while not IsOne(xx) do
      xx := xx^a[1];
      k := k-1;
    od;
    res := res/a[1]^k;
  od;
  return res;
end);



