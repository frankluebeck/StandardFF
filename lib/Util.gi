#############################################################################
##
#A  Util.gi                                                      Frank Lübeck
##  
##  The files Util.g{i,d}  contain code for some utiltiy functions.
##  


##  <#GAPDoc Label="SimpleRandomRange">
##  <ManSection>
##  <Func Name="SimpleRandomRange" Arg="max, seed" />
##  <Returns>a pair <C>[rand, nseed]</C> of integers </Returns>
##  <Description>
##  The argument <A>max</A> must be  a non-negative integer and <A>seed</A> an
##  integer  which is  no  multiple  of the  prime  <M>2^{32}-5</M> (only  its
##  residue class modulo this prime is  used). This function implements a very
##  simple pseudo  random number generator defined  in <Cite Key="StdFFCyc"/>.
##  The result  is a pair <C>[rand,  nseed]</C> where <C>rand</C> is  a random
##  integer  in  the range  <C>[0..<A>max</A>-1]</C>  and  <C>nseed</C> is  an
##  integer that can be used as seed for further calls of this function.
##  <P/> 
##  <Example>gap> pair := [0,1];;
##  gap> for i in [1..3] do
##  >   pair := SimpleRandomRange(1000000, pair[2]); Print(pair, "\n");
##  > od;
##  [ 313679, 1347244577 ]
##  [ 625568, 2686795180 ]
##  [ 853364, 3665171566 ]
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# seed must be <> 0 mod 2^32-5.
# max some positive integer
# returns:  [rand, newseed] with 0 <= rand < max and new seed
PROPERRANDOMRANGE := false;
InstallGlobalFunction(SimpleRandomRange, function(max, seed)
  local m, a, b, r, mm;
  if PROPERRANDOMRANGE = true then
    return [Random(0,max-1),1];
  fi;
  # prime 2^32-5;
  m := 4294967291;
  # has order m-1 mod m
  a := 1347244577;
  # compute random (m-1)-adic digits with maximum >= 100*max
  b := max*100;
  # pseudo random number r and possible maximum mm
  r := 0;
  mm := 1;
  while mm < b do
    mm := mm*(m-1);
    seed := a*seed mod m;
    r := r*(m-1) + seed - 1;
  od;
  return [QuoInt(max*r, mm), seed];
end);

##  <#GAPDoc Label="FindLinearCombination">
##  <ManSection>
##  <Func Name="FindLinearCombination" Arg="v, start"/>
##  <Returns> a pair <C>[serec, lk]</C> of a record and vector 
##  or <K>fail</K></Returns>
##  <Description>
##  Repeated calls  of this  function build  up a  semiechelon basis  from the
##  given  arguments <A>v</A>  which  must  be row  vectors.  To initialize  a
##  computation  the function  is  called  with a  start  vector <A>v</A>  and
##  <K>false</K> as  second argument.  The return value  is a  pair <C>[serec,
##  lk]</C>  where <C>serec</C>  is  a  record which  collects  data from  the
##  previous  calls of  the  function  and <C>lk</C>  is  a  row vector  which
##  expresses  <A>v</A> as  linear combination  of the  vectors from  previous
##  calls,  or <K>fail</K>  if there  is no  such linear  combination. In  the
##  latter  case  the  data  in  the record  is  extended  with  the  linearly
##  independent vector <C>v</C>.
##  <P/>
##  In the following example  we show how to compute a  divisor of the minimal
##  polynomial of a matrix.
##  <Example>gap> mat := Product(GeneratorsOfGroup(Sp(30,5)));;
##  gap> x := Indeterminate(GF(5), "x");;
##  gap> v := (mat^0)[1];;
##  gap> b := FindLinearCombination(v, false);;
##  gap> repeat
##  >   v := v*mat;
##  >   l := FindLinearCombination(v, b[1]);
##  > until IsList(l[2]);
##  gap> mp := Value(UnivariatePolynomial(GF(5),
##  >         Concatenation(-l[2], [One(GF(5))])), x);
##  x^30+Z(5)^3*x^29+Z(5)^3*x+Z(5)^0
##  gap> # equal to minimal polynomial because of degree
##  gap> mp = Value(MinimalPolynomial(GF(5), mat), x);
##  true
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallGlobalFunction(FindLinearCombination, function(v, start)
  local lc, c, pos, i;
  if start = false then
    start := rec(orig := [],
                 sebasis := [],
                 lincomb := [],
                 pivots := []);
  fi;
  if ForAll(v, IsZero) then
    if Length(start.orig) = 0 then
      return [];
    else
      return List(start.orig, a-> 0*a[1]);
    fi;
  fi;
  # reduce
  Add(start.orig, v);
  lc := List(start.sebasis, a-> 0*v[1]);
  Add(lc, One(v[1]));
  for i in [1..Length(start.sebasis)] do
    c := v[start.pivots[i]];
    if not IsZero(c) then
      v := v-c*start.sebasis[i];
      lc := lc-c*start.lincomb[i];
    fi;
  od;
  if ForAll(v, IsZero) then
    Remove(start.orig);
    return [start, -lc{[1..Length(start.sebasis)]}];
  fi;
  pos := First([1..Length(v)], i-> not IsZero(v[i]));
  c := v[pos]^-1;
  v := c*v;
  lc := c*lc;
  Add(start.sebasis, v);
  Add(start.lincomb, lc);
  Add(start.pivots, pos);
  return [start, fail];
end);

InstallMethod(MinimalPolynomial, 
                 ["IsField and IsFinite", "IsAlgebraicElement", "IsPosInt"],
function(F, a, nr)
  local fam, d, f, z, mat, c, lc, m, i;
  fam := FamilyObj(a);
  if not F = fam!.baseField then
    TryNextMethod();
  fi;
  d := fam!.deg;
  f := fam!.polCoeffs;
  f := - f[d+1]^-1*f{[1..d]};
  z := Zero(F);
  # mat becomes the matrix of a with respect to the standard
  # basis 1, X, X^2, ..., X^{d-1} of the algebraic extension
  mat := NullMat(d,d,F);
  c := ShallowCopy(ExtRepOfObj(a));
  lc := Length(c);
  mat[1]{[1..lc]} := c;
  for i in [2..d] do
    if lc=d then
      m := c[d];
      c{[2..d]} := c{[1..d-1]};
      c[1] := z;
      if m <> z then
        c := c + m*f;
      fi;
    else
      c{[2..lc+1]} := c{[1..lc-1]};
      c[1] := z;
      lc := lc+1;
    fi;
    mat[i]{[1..lc]} := c;
  od;
  return MinimalPolynomial(F, mat, nr);
end);

InstallOtherMethod(MinimalPolynomial, 
                 ["IsField and IsFinite", "IsAlgebraicElement",
                 "IsUnivariatePolynomial"],
function(F, a, var)
  return MinimalPolynomial(F, a,
                    IndeterminateNumberOfUnivariateRationalFunction(var));
end);

##  <#GAPDoc Label="BerlekampMassey">
##  <ManSection>
##  <Func Name="BerlekampMassey" Arg="u" />
##  <Returns>a list of field elements</Returns>
##  <Description>
##  The argument <A>u</A> is a list of elements in a field <M>F</M>.
##  This function implements the Berlekamp-Massey algorithm which returns
##  the shortest sequence <M>c</M> of elements in <M>F</M> such that for each
##  <M>i > l</M>, the length of <M>c</M>, we have <M>u[i] = \sum_{{j=1}}^l 
##  <A>u</A>[i-j] c[j]</M>.
##  <Example>gap> x := Indeterminate(GF(23), "x");;
##  gap> f := x^5 + Z(23)^16*x + Z(23)^12;;
##  gap> u := List([1..50], i-> Value(x^i mod f, 0));;
##  gap> c := BerlekampMassey(u);;
##  gap> ForAll([6..50], i-> u[i] = Sum([1..5], j-> u[i-j]*c[j]));
##  true
##  gap> -c;
##  [ 0*Z(23), 0*Z(23), 0*Z(23), Z(23)^16, Z(23)^12 ]
##  </Example>
##  </Description>
##  </ManSection>
##  
##  <ManSection>
##  <Func Name="MinimalPolynomialByBerlekampMassey" Arg="x" />
##  <Returns>the minimal polynomial of <A>x</A></Returns>
##  <Description>
##  Here <M>x</M> must be an element of an algebraic extension field <M>F/K</M>.
##  (<M>K</M> must be the <Ref BookName="Reference" Attr="LeftActingDomain" />
##  of <M>F</M>). 
##  This function computes the minimal polynomial of <A>x</A> over <M>K</M> by
##  applying the Berlekamp-Massey algorithm to the list of traces of <M><A>x</A>^i</M>.
##  <Example>gap> x := Indeterminate(GF(23), "x");;
##  gap> f := x^5 + Z(23)^16*x + Z(23)^12;;
##  gap> F := AlgebraicExtension(GF(23), f);;
##  gap> mp := MinimalPolynomialByBerlekampMassey(PrimitiveElement(F));;
##  gap> Value(mp, x) = f;
##  true
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
InstallGlobalFunction(BerlekampMassey, function(u)
  local N, l, a, r, t, b, z, dd, d, c, a1, m, n, k;
  N := Length(u);
  l := 0;
  a := [];
  r := -1;
  t := 0;
  b := [];
  z := Zero(u[1]);
  dd := One(u[1]);
  for n in [0..N-1] do
    d := u[n+1];
    for k in [1..l] do
      d := d-a[k]*u[n-k+1];
    od;
    if not IsZero(d) then
      c := d*dd^-1;
      a1 := a - c*Concatenation(z*[1..n-r],b);
      a1[n-r] := a1[n-r] + c;
      if 2*l <= n then
        m := n+1-l;
        t := l;
        l := m;
        b := a;
        r := n;
        dd := d;
      fi;
      a := a1;
    fi;
  od;
  return a;
end);

InstallGlobalFunction(MinimalPolynomialByBerlekampMassey,  function(x)
  local cl, F, K, n, y, l, res, i;
  cl := function(y)
    res := y![1];
    if not IsList(res) then
      res := [res];
    fi;
    return res;
  end;
  F := DefaultField(x);
  K := LeftActingDomain(F);
  n := Degree(DefiningPolynomial(F));
  y := x^0;
  l := [cl(y)[1]];
  for i in [1..2*n-1] do
    y := y*x;
    Add(l, cl(y)[1]);

  od;
  res := -Reversed(BerlekampMassey(l));
  ConvertToVectorRep(res);
  Add(res, One(res[1]));
  return UnivariatePolynomial(K,res,1);
end);

# with multiplicative upper bound
# (this is sufficient in this context, 'OrderMod' would try lengthy
# or impossible integer factorizations for larger parameters)
InstallGlobalFunction(OrderModBound, function(a, m, bound) 
  local d;
  d := DivisorsInt(bound);
  return First(d, n-> PowerMod(a, n, m) = 1);
end);

##  <#GAPDoc Label="DLog">
##  <ManSection>
##  <Func Name="DLog" Arg="base, x, m"/>
##  <Returns>an integer</Returns>
##  <Description>
##  The argument  <A>base</A> must  be a  multiplicative element  and <A>x</A>
##  must lie in the cyclic group  generated by <A>base</A>. The third argument
##  <A>m</A>  must be  the order  of  <A>base</A> or  its factorization.  This
##  function returns the discrete logarithm,  that is an integer <M>e</M> such
##  that <A>base</A><M>^e = </M> <A>x</A>.
##  <P/> 
##  If  <A>m</A>  is  prime  then  Shanks'  algorithm  is  used  (which  needs
##  <M>O(\sqrt{<A>m</A>})</M> space and time). Otherwise  let <A>m</A> <M> = r
##  l</M> and  <M>e =  a +  b r</M>  with <M>0  \leq a  &lt; r</M>.  Then <M>a
##  =</M>  <C>DLog</C><M>(<A>base</A>^l, <A>x</A>^l,  r)</M> and  <M>b =  </M>
##  <C>DLog</C><M>(<A>base</A>^r, <A>x</A>/<A>base</A>^a, l)</M>.
##  <P/>
##  This  function  is  used  for   a  method    of  <Ref  BookName="Reference"
##  Oper="LogFFE"/>.
##  
##  <Example>gap> F := FF(67, 12);
##  FF(67, 12)
##  gap> st := SteinitzPairConwayGenerator(F);
##  [ 12, 5575571448927120404890 ]
##  gap> z := ElementSteinitzNumber(F, st[2]);;
##  gap> x := StandardPrimitiveRoot(F);;
##  gap> DLog(z, x, Size(F)-1);
##  103390307992785367583
##  gap> K := GF(67,12);
##  GF(67^12)
##  gap> zz := Z(67^12);
##  z
##  gap> LogFFE(zz^2+1, zz);
##  1667375214152688471247
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# helper: returns e <= r such that base^e = x or fail
InstallGlobalFunction(DLogShanks, function(base, x, r)
  local rr, baby, ord, giant, t, pos, i, j;
  rr := RootInt(r, 2);
  baby := [One(base)];
  if x = baby[1] then
    return 0;
  fi;
  for i in [1..rr-1] do
    baby[i+1] := baby[i]*base;
    if x = baby[i+1] then
      return i;
    fi;
  od;
  giant := baby[rr]*base;
  ord := [0..rr-1];
  SortParallel(baby, ord);
  t := x;
  for j in [1..QuoInt(r, rr)+1] do
    t := t*giant;
    pos := PositionSet(baby, t);
    if IsInt(pos) then
      return (ord[pos] - j * rr) mod r;
    fi;
  od;
  return fail;
end);
# recursive method, m can be order of base or its factorization
# Let r be the largest prime factor of m, then we use
#     base^e = x with e = a + b*r where 0 <= a < r and
# 0 <= b < m/r, and compute a with DLogShanks and b by
# recursion.
InstallGlobalFunction(DLog, function(base, x, m)
  local r, mm, mp, a, b;
  if not IsList(m) then
    m := Factors(m);
  fi;
  if Length(m) = 1 then
    return DLogShanks(base, x, m[1]);
  fi;
  r := m[Length(m)];
  mm := m{[1..Length(m)-1]};
  mp := Product(mm);
  a := DLogShanks(base^mp, x^mp, r);
  b := DLog(base^r, x/(base^a), mm);
  return a + b*r;
end);

# this seems to be better than the GAP library function in many cases.
InstallMethod(LogFFE, ["IsFFE and IsCoeffsModConwayPolRep", 
                       "IsFFE and IsCoeffsModConwayPolRep"], 
function(x, base)
  local ob, o, e;
  ob := Order(base);
  o := Order(x);
  if  ob mod o <> 0 then
    return fail;
  fi;
  if ob <> o then
    e := ob/o;
    base := base^e;
  else
    e := 1;
  fi;
  return DLog(base, x, o) * e;
end);

##  <#GAPDoc Label="FindConjugateZeroes">
##  <ManSection>
##  <Func Name="FindConjugateZeroes" Arg="K, cpol, qq"/>
##  <Returns>a list of field elements </Returns>
##  <Description>
##  The arguments  must be a  finite field <A>K</A>, a  polynomial <A>cpol</A>
##  over  <A>K</A>  (or its  coefficient  list)  and  the order  <A>qq</A>  of
##  a  subfield  of   <A>K</A>.  The  polynomial  must   have  coeffcients  in
##  the  subfield  with <A>qq</A>  elements,  must  be irreducible  over  this
##  subfield and  split into linear  factors over <A>K</A>. The  function <Ref
##  Func="FindConjugateZeroes"/> returns the list  of zeroes of <A>cpol</A> in
##  <A>K</A>.
##  <Example>gap> K := GF(67,18);
##  GF(67^18)
##  gap> F := FF(67,18);
##  FF(67, 18)
##  gap> p1 := DefiningPolynomial(K);;
##  gap> p2 := DefiningPolynomial(F);;
##  gap> lK := FindConjugateZeroes(K, p2, 67);;
##  gap> lF := FindConjugateZeroes(F, p1, 67);;
##  gap> Minimum(List(lF, SteinitzNumber));
##  12578642085854203709864918398540
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# K finite field, cpol polynomial over K irreducible over subfield with qq
# elements and zeroes in K.
InstallGlobalFunction(FindConjugateZeroes, function(K, cpol, qq)
  local fac, d, r, rr, g1, g2, z, i;
  if IsPolynomial(cpol) then
    cpol := One(K)*CoefficientsOfUnivariatePolynomial(cpol);
  fi;
  if Characteristic(K) = 2 then
    return FindConjugateZeroesChar2(K, cpol, qq);
  fi;
  fac := cpol;
  d := Length(fac)-1;
  while Length(fac) > 2 do
    r := 0;
    while IsZero(r) do
      r := [Random(K),One(K)];
    od;
    rr := PowerModCoeffs(r, (Size(K)-1)/2, fac);
    rr[1] := rr[1]-One(K);
    g1 := GcdCoeffs(rr, fac);
    rr[1] := rr[1] + 2*One(K);
    g2 := GcdCoeffs(rr, fac);
    if Length(g1) > Length(g2) and Length(g2) > 1 then
      fac := g2;
    elif Length(g1) > 1 then
      fac := g1;
    fi;
  od;
  z := [-fac[1]/fac[2]];
  for i in [1..d-1] do
    Add(z, z[i]^qq);
  od;
  return z;
end);

##  <#GAPDoc Label="ZeroesConway">
##  <ManSection>
##  <Func Name="ZeroesConway" Arg="F" />
##  <Returns>a list of field elements </Returns>
##  <Description>
##  Here, <A>F</A>  must be a  standard finite  field, say of  degree <M>n</M>
##  over the  prime field  with <M>p</M> elements.  This function  returns the
##  same as <C>FindConjugateZeroes(F, One(F)*ConwayPol(p,  n), p)</C> (using a
##  specific implementation).
##  <P/> 
##  <Example>gap> F := FF(23,29);
##  FF(23, 29)
##  gap> l := Set(FindConjugateZeroes(F, One(F)*ConwayPol(23,29), 23));;
##  gap> l = Set(ZeroesConway(F));
##  true
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# This is a specialized version of FindConjugateZeroes for the Conway
# polynomial of the field K. It returns the list of all zeroes in K.
InstallGlobalFunction(ZeroesConway, function(K)
  local p, n, f, x, l, bas, fac, len, pr, c, rr, g1, g2, z, i, a;
  p := Characteristic(K);
  n := DegreeOverPrimeField(K);
  # Conway polynomial as polynomial coeffs over  K
  fac := One(K)*ConwayPol(p, n);
  if p = 2 then
    return FindConjugateZeroesChar2(K, fac, 2);
  fi;
  len := Length(fac);
  f := GF(p^n);
  # we precompute for y = x^((p-1)/2) the list of y^(p^i) mod Conway polynomial,
  # i < n (using GAPs arithemetic in GF(p^n)).
  x := PrimitiveElement(f);
  l := [x^((p-1)/2)];
  for i in [1..n-1] do
    Add(l, l[i]^p);
  od;
  if IsPositionalObjectRep(l[n]) then
    l := List(l, a->a![1]);
  else
    bas := CanonicalBasis(f);
    l := List(l, a-> One(K)*Coefficients(bas, a));
  fi;
  l := One(K)*l;
  pr := function(l1, l2)
    local ll;
    ll := ProductCoeffs(l1,l2);
    ReduceCoeffs(ll, fac);
    ShrinkRowVector(ll);
    return ll;
  end;
  # now we split fac
  while len > 2 do
    # This is Cantor-Zassenhaus, since the factors of fac over K are linear 
    # it is sufficient to try random polynomials of form X+c, c in K.
    # We compute (X + c) ^ ((p^n-1)/2) mod fac using
    #      (X+c)^(p^i) = X^(p^i) + c^(p^i)
    # (first summands mod fac computed above, second is computation in K)
    # and
    #      (p^n-1)/2 = (p-1)/2 * (1 + p + p^2 + ...+ p^(n-1))
    c := Random(K);
    rr := ShallowCopy(l[1]);
    rr[1] := rr[1] + c;
    for i in [2..n] do
      c := c^p;
      l[i][1] := l[i][1] + c;
      rr := pr(rr, l[i]);
      l[i][1] := l[i][1] - c;
    od;
    # now Gcd of this power +/- 1 with fac is likely to split fac
    rr[1] := rr[1]-One(K);
    g1 := GcdCoeffs(rr, fac);
    rr[1] := rr[1] + 2*One(K);
    g2 := GcdCoeffs(rr, fac);
    if Length(g1) > Length(g2) and Length(g2) > 1 then
      fac := g2;
    elif Length(g1) > 1 then
      fac := g1;
    fi;
    if Length(fac) < len then
      # in case of splitting we can from now compute modulo found factor fac
      for a in l do
        ReduceCoeffs(a, fac);
        ShrinkRowVector(a);
      od;
      len := Length(fac);
    fi;
  od;
  z := [-fac[1]/fac[2]];
  for i in [1..n-1] do
    Add(z, z[i]^p);
  od;
  return z;
end);

InstallGlobalFunction(FindConjugateZeroesChar2, function(K, cpol, qq)
  local fac, d, m, r, rr, g1, g2, z, j, i;
  fac := cpol;
  d := Length(fac)-1;
  m := DegreeOverPrimeField(K);
  while Length(fac) > 2 do
    r := [Random(K), Random(K)];
    rr := r;
    for j in [2..m] do
      r := PowerModCoeffs(r, 2, fac);
      rr := rr + r;
    od;
    g1 := GcdCoeffs(rr, fac);
    rr[1] := rr[1] + One(K);
    g2 := GcdCoeffs(rr, fac);
    if Length(g1) > Length(g2) and Length(g2) > 1 then
      fac := g2;
    elif Length(g1) > 1 then
      fac := g1;
    fi;
  od;
  z := [-fac[1]/fac[2]];
  for i in [1..d-1] do
    Add(z, z[i]^qq);
  od;
  return z;
end);

# extend list of Brent factors in FactInt:

# this is essentially FetchBrentFactors from FactInt
# application examples (only needed once per GAP installation):
### Brent original factor lists, no longer maintained
##  FetchMoreFactors(
##    "http://wwwmaths.anu.edu.au/~brent/ftp/factors/factors.gz",
##    false);
### newer and more factors collected by Jonathan Crombie, check occasionally
### if newer updates are available and adjust date part of URL:
##  FetchMoreFactors(
##    "http://myfactorcollection.mooo.com:8090/brentdata/Feb28_2021/factors.gz",
##    true);
InstallGlobalFunction(FetchMoreFactors,   function ( url, write )
  local  str, get, comm, rows, b, k, a, dir;

  # Fetch the file from R. P. Brent's ftp site and gunzip it into 'str'.

  str := "";
  get := OutputTextString(str, false);
  comm := Concatenation("wget -q ", url, " -O - | gzip -dc ");
  Process(DirectoryCurrent(), Filename(DirectoriesSystemPrograms(),"sh"),
          InputTextUser(), get, ["-c", comm]);

  rows := SplitString(str, "", "\n");
  str := 0;
  for a in rows do
    b := List(SplitString(a, "", "+- \n"), Int);
    if not IsBound(BRENTFACTORS[b[1]]) then
      BRENTFACTORS[b[1]] := [];
    fi;
    if '-' in a then
      k := b[2];
    else
      k := 2*b[2];
    fi;
    if not IsBound(BRENTFACTORS[b[1]][k]) then
      BRENTFACTORS[b[1]][k] := [b[3]];
    else
      Add(BRENTFACTORS[b[1]][k], b[3]);
    fi;
  od;
  if write then
    dir := GAPInfo.PackagesInfo.("factint")[1].InstallationPath;
    WriteBrentFactorsFiles(Concatenation(dir,"/tables/brent/"));
  fi;
end);

##  <#GAPDoc Label="StandardValuesBrauerCharacter">
##  <ManSection>
##  <Func Name="StandardValuesBrauerCharacter" Arg="tab, bch"/>
##  <Returns>a Brauer character </Returns>
##  <Func Name="IsGaloisInvariant" Arg="tab, bch"/>
##  <Returns><K>true</K> or <K>false</K></Returns>
##  <Description>
##  The argument  <A>tab</A> must be  a Brauer  character table for  which the
##  Brauer characters  are defined with  respect to  the lift given  by Conway
##  polynomials. And  <A>bch</A> must  be an  irreducible Brauer  character of
##  this table.
##  <P/> 
##  The   function   <Ref  Func="StandardValuesBrauerCharacter"/>   recomputes
##  the    values    corresponding    to    the    lift    given    by    <Ref
##  Func="StandardCyclicGenerator"/>, provided that the Conway polynomials for
##  computing the Frobenius  character values of <A>bch</A>  are available. If
##  Conway  polynomials are  missing  the corresponding  character values  are
##  substituted by <K>fail</K>. If the  result does not contain <K>fail</K> it
##  is a  class function  which is  Galois conjugate  to <A>bch</A>  (see <Ref
##  BookName="Reference" Meth="GaloisCyc" Label="for a class function"/>).
##  <P/>
##  The  utility <Ref  Func="IsGaloisInvariant"/> returns  <K>true</K> if  all
##  Galois conjugates  of <A>bch</A> are  Brauer characters in  <A>tab</A>. If
##  this is the  case then different lifts will permute  the Galois conjugates
##  and all of them are Brauer characters with respect to any lift.
##  
##  <Example>gap> tab := BrauerTable("M", 19);
##  BrauerTable( "M", 19 )
##  gap> # cannot translate some values to different lift
##  gap> fail in StandardValuesBrauerCharacter(tab, Irr(tab)[16]);
##  false
##  gap> # but table contains the irreducible Brauer characters for any lift
##  gap> ForAll(Irr(tab), bch-> IsGaloisInvariant(tab, bch));
##  true
##  gap> tab := BrauerTable("A18", 7);
##  BrauerTable( "A18", 7 )
##  gap> # here different lifts lead to different Brauer character tables
##  gap> bch := Irr(tab)[123];;
##  gap> IsGaloisInvariant(tab, bch);
##  false
##  gap> new := StandardValuesBrauerCharacter(tab, bch);;
##  gap> fail in new;
##  false
##  gap> Position(Irr(tab), new);
##  fail
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# translate Brauer character with respect to Conway lifting to standard lifting
InstallGlobalFunction(StandardValuesBrauerCharacter, function(tab, bch)
  local o, p, res, ev, ord, k, f, base, sp, z, c, e, l, i, j;
  o := OrdersClassRepresentatives(tab);
  p := UnderlyingCharacteristic(tab);
  res := [];
  for i in [1..Length(o)] do
    if 2 mod o[i] = 0 or Conductor(bch[i]) = 1 then
      res[i] := bch[i];
    else
      ev := EigenvaluesChar(tab, bch, i);
      ord := Length(ev);
      k := OrderMod(p, ord);
      f := FF(p,k);
      base := StandardCyclicGenerator(f, ord);
      sp := SteinitzPairConwayGenerator(f);
      if sp = fail then
        res[i] := fail;
      else
        z := ElementSteinitzNumber(f, sp[2]);
        c := (p^k-1)/ord;
        e := DLog(base, z^c, ord);
        l := [E(ord)^e];
        for j in [1..ord-1] do
          Add(l, l[j]*l[1]);
        od;
        res[i] := ev*l;
      fi;
    fi;
  od;
  return Character(tab, res);
end);
InstallGlobalFunction(IsGaloisInvariant, function(tab, bch)
  local N;
  N := Conductor(bch);
  return ForAll(PrimeResidues(N), j-> GaloisCyc(bch, j) in Irr(tab));
end);


