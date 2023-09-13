#############################################################################
##
#A  Util.gi                                                      Frank Lübeck
##  
##  The files Util.g{i,d}  contain code for some utiltiy functions.
##  


##  <#GAPDoc Label="InvModCoeffs">
##  <ManSection>
##  <Oper Name="InvModCoeffs" Arg="fcoeffs, gcoeffs" />
##  <Returns>a list of <K>fail</K></Returns>
##  <Description>
##  The arguments <A>fcoeffs</A> and <A>gcoeffs</A> are coeffient lists of two
##  polynomials <M>f</M> and <M>g</M>.
##  This operation returns  the coefficient list of  the inverse <M>f^{-1}</M>
##  modulo  <M>g</M>, if  <M>f</M> and  <M>g</M> are  coprime, and <K>fail</K>
##  otherwise.
##  <P/>
##  The  default  method  computes  the  inverse  by  the  extended  Euclidean
##  algorithm.
##  <Example>gap> f := Z(13)^0*[ 1, 10, 1, 11, 0, 1 ];;
##  gap> g := Z(13)^0*[ 5, 12, 5, 12, 2, 0, 2 ];;
##  gap> InvModCoeffs(f, g);
##  fail
##  gap> GcdCoeffs(f, g);
##  [ Z(13)^0, 0*Z(13), Z(13)^0 ]
##  gap> f[1]:=f[1]+1;;
##  gap> finv := InvModCoeffs(f, g);
##  [ Z(13)^9, Z(13)^10, Z(13)^10, Z(13)^8, Z(13)^5, Z(13)^6 ]
##  gap> pr := ProductCoeffs(finv, f);;
##  gap> ReduceCoeffs(pr, g);; ShrinkRowVector(pr);; pr;
##  [ Z(13)^0 ]
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
# f^-1 mod g
InstallMethod(InvModCoeffs, [IsList, IsList], function(f, g)
  local rf, z, rg, t, h, rh, i;
  rf:=[One(g[1])];
  z:=Zero(g[1]);
  rg:=[];
  while g<>[] do
    t:=QuotRemPolList(f,g);
    h:=g;
    rh:=rg;
    g:=t[2];
    if Length(t[1])=0 then
      rg:=[];
    else
      rg:=ShallowCopy(-ProductCoeffs(t[1],rg));
    fi;
    for i in [1..Length(rf)] do
      if IsBound(rg[i]) then
        rg[i]:=rg[i]+rf[i];
      else
        rg[i]:=rf[i];
      fi;
    od;
    f:=h;
    rf:=rh;
    ShrinkRowVector(g);
  od;
  if Length(f) <> 1 then
    return fail;
  else
    return 1/f[1]*rf;
  fi;
end);

##  <#GAPDoc Label="StandardAffineShift">
##  <ManSection>
##  <Func Arg="q, i" Name="StandardAffineShift" />
##  <Returns>an integer in range <C>[0..q-1]</C></Returns>
##  <Description>
##  This function returns <M>(m <A>i</A>  + a)  \textrm{ mod } <A>q</A></M>,
##  where <M>m</M>  is the largest integer  prime to <A>q</A> and  <M>\leq 4
##  <A>q</A> /  5</M>, and  a is  the largest integer  <M>\leq 2  <A>q</A> /
##  3</M>.
##  <P/>
##  For  fixed <M>q</M>  this function  provides  a bijection  on the  range
##  <C>[0..q-1]</C>.
##  <Example>gap> List([0..10], i-> StandardAffineShift(11, i));
##  [ 7, 4, 1, 9, 6, 3, 0, 8, 5, 2, 10 ]
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
BindGlobal("SASCache", [0,0,0]);
InstallGlobalFunction(StandardAffineShift, function(q, i)
  local m, a;
  if q = SASCache[1] then
    m := SASCache[2];
    a := SASCache[3];
  else 
    m := QuoInt(4*q, 5);
    while Gcd(m, q) <> 1 do
      m := m-1;
    od;
    a := QuoInt(2*q, 3);
    SASCache{[1..3]} := [q,m,a];
  fi;
  return (m*i+a) mod q;
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
##  The argument  <A>u</A> is  a list  of elements in  a field  <M>F</M>. This
##  function  implements  the  Berlekamp-Massey algorithm  which  returns  the
##  shortest sequence <M>c</M> of elements in <M>F</M> such that for each <M>i
##  >  l</M>,  the  length  of  <M>c</M>, we  have  <M>u[i]  =  \sum_{{j=1}}^l
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
##  <Func Name="MinimalPolynomialByBerlekampMasseyShoup" Arg="x" />
##  <Returns>the minimal polynomial of <A>x</A></Returns>
##  <Description>
##  Here   <M>x</M>   must  be   an   element   of  an   algebraic   extension
##  field  <M>F/K</M>.   (<M>K</M>  must  be  the   <Ref  BookName="Reference"
##  Attr="LeftActingDomain"  />  of  <M>F</M>).  This  function  computes  the
##  minimal   polynomial   of  <A>x</A>   over   <M>K</M>   by  applying   the
##  Berlekamp-Massey algorithm to the list of traces of <M><A>x</A>^i</M>.
##  <P/>
##  The second variant uses the algorithm by Shoup in <Cite Key="ShoupMiPo"/>.
##  <Example>gap> x := Indeterminate(GF(23), "x");;
##  gap> f := x^5 + Z(23)^16*x + Z(23)^12;;
##  gap> F := AlgebraicExtension(GF(23), f);;
##  gap> mp := MinimalPolynomialByBerlekampMassey(PrimitiveElement(F));;
##  gap> Value(mp, x) = f;
##  true
##  gap> mp = MinimalPolynomialByBerlekampMasseyShoup(PrimitiveElement(F));
##  true
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
InstallGlobalFunction(BerlekampMassey, function(u)
  local N, l, a, r, t, b, z, ur, dd, d, c, a1, m, n;
  N := Length(u);
  l := 0;
  a := [];
  r := -1;
  t := 0;
  b := [];
  z := 0*u;
  ur := Reversed(u);
  dd := One(u[1]);
  for n in [0..N-1] do
    d := u[n+1];
##      for k in [1..l] do
##        d := d-a[k]*u[n-k+1];
##      od;
    if l > 0 then
      d := d - a * ur{[N+1-n..N+l-n]};
    fi;
    if not IsZero(d) then
      c := d*dd^-1;
      a1 := a - Concatenation(z{[1..n-r]},c*b);
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

# more efficient and more complicated version by Shoup
# first a utility function for transposed multiplication
# tau is field element, v vector, returns tau \circ v
# see [Shoup, Efficient Computation of Minimal Polynomials, Sec. 3]
BindGlobal("TransposedMult", function(tau, v)
  local fam, f, n, b, ft, xx, h, bt, t1, t2, t3, res, zv;
  fam := FamilyObj(tau);
  f := fam!.polCoeffs;
  n := Length(f)-1;
  b := tau![1];
  if not IsList(b) then
    b := [b];
  fi;
  while Length(b) < n do
    Add(b, Zero(fam));
  od;
  bt := Reversed(b);
  zv := Zero(f{[1]});
  ft := Reversed(f);
  if IsBound(fam!.ftinv) then
    h := fam!.ftinv;
  else
    xx := 0*ShallowCopy(ft);
    Remove(xx);
    xx[n] := One(f[1]);
    h := InvModCoeffs(ft, xx);
    ConvertToVectorRep(h);
    fam!.ftinv := h;
  fi;
  t1 := ProductCoeffs(ft, v);
  t1 := t1{[n+1..Length(t1)]};
  t2 := ProductCoeffs(bt, v);
  t2 := t2{[n..Length(t2)]};
  t3 := ProductCoeffs(h,bt);
  if Length(t3) > n-1 then
    t3 := t3{[1..n-1]};
  fi;
  t3 := ProductCoeffs(t3, t1);
  if Length(t3) > n-1 then
    t3 := t3{[1..n-1]};
  fi;
  t3 := Concatenation(zv, t3);
  res := t2 - t3;
  ConvertToVectorRep(res);
  return res;
end);

# this is like TransposedMult(tau, v) but with some precomputed
# data (good, when called multiple times with the same tau):
#    ft : coeffs of reversed (deg n) defining polynomial of extension
#         containig tau
#    bt : coeffs of reversed (deg n-1) of tau
#    hbt : coeffs of h * bt mod X^{n-1}, where h is inverse of
#          ft mod X^{n-1}
BindGlobal("TransposedMultCoeffsPre", function(v, ft, bt, hbt)
  local n, zv, t1, t2, t3, res;
  n := Length(ft)-1;
  zv := Zero(ft{[1]});
  t1 := ProductCoeffs(ft, v);
  t1 := t1{[n+1..Length(t1)]};
  t2 := ProductCoeffs(bt, v);
  t2 := t2{[n..Length(t2)]};
  t3 := ProductCoeffs(hbt, t1);
  if Length(t3) > n-1 then
    t3 := t3{[1..n-1]};
  fi;
  t3 := Concatenation(zv, t3);
  res := t2 - t3;
  ConvertToVectorRep(res);
  return res;
end);

# simple version for field extensions of
# [Shoup, Efficient Computation of Minimal Polynomials, Sec. 4]
InstallGlobalFunction(MinimalPolynomialByBerlekampMasseyShoup, function(x)
  local cl, F, K, n, k, yy, z, bt, ft, fam, h, xx, hbt, v, l, res, i, j;
  # field element as coefficient list
  cl := function(y)
    local res;
    res := y![1];
    if not IsList(res) then
      res := [res];
    fi;
    return res;
  end;
  F := DefaultField(x);
  K := LeftActingDomain(F);
  n := Degree(DefiningPolynomial(F));
  if n=1 then
    res := [-x, One(K)];
    return UnivariatePolynomial(K,res,1);
  fi;

  # to be tuned: number of precomputed powers of x
  k := 2*RootInt(n,2);
  yy := [x^0, x];
  for i in [2..k] do
    Add(yy, yy[i]*x);
  od;
  for i in [1..k+1] do
    yy[i] := cl(yy[i]);
  od;

  # x^k is needed for transposed multiplication
  # precomputation:
  z := Zero(K);
  while Length(yy[k+1]) < n do
    Add(yy[k+1], z);
  od;
  bt := Reversed(yy[k+1]);
  ft := Reversed(CoefficientsOfUnivariatePolynomial(DefiningPolynomial(F)));
  # we cache inverse of ft modulp X^{n-1}
  fam := FamilyObj(x);
  if IsBound(fam!.ftinv) then
    h := fam!.ftinv;
  else
    xx := 0*ShallowCopy(ft);
    Remove(xx);
    xx[n] := One(K);
    h := InvModCoeffs(ft, xx);
    ConvertToVectorRep(h);
    fam!.ftinv := h;
  fi;
  # h * bt mod X^{n-1}
  hbt := ProductCoeffs(h,bt);
  if Length(hbt) > n-1 then
    hbt := hbt{[1..n-1]};
  fi;

  # now we can compute efficiently the absolute coeffs of
  # x^i for i = 0 .. 2n-1
  v := ShallowCopy(cl(x){[1]});
  v[1] := One(K);
  l := [];
  while Length(l) < 2*n do
    for j in [1..k] do
      Add(l, yy[j]*v);
    od;
    if Length(l) < 2*n then
      v := TransposedMultCoeffsPre(v, ft, bt, hbt);
    fi;
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

if not IsBound(DLog) then
  # otherwise this code was moved to and is already read from the GAP library

##  <#GAPDoc Label="DLog">
##  <ManSection>
##  <Func Name="DLog" Arg="base, x[, m]"/>
##  <Returns>an integer</Returns>
##  <Description>
##  The argument  <A>base</A> must  be a  multiplicative element  and <A>x</A>
##  must lie in the cyclic group  generated by <A>base</A>. The third argument
##  <A>m</A>  must  be the  order  of  <A>base</A>  or its  factorization.  If
##  <A>m</A> is  not given, it  is computed  first. This function  returns the
##  discrete logarithm, that is an integer <M>e</M> such that <A>base</A><M>^e
##  = </M> <A>x</A>.
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
##  [ 12, 5118698034368952035290 ]
##  gap> z := ElementSteinitzNumber(F, st[2]);;
##  gap> x := StandardPrimitiveRoot(F);;
##  gap> DLog(z, x, Size(F)-1);
##  231901568073107448223
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
BindGlobal("DLogShanks", function(base, x, r)
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
BindGlobal("DLog", function(base, x, m...)
  local r, mm, mp, a, b;
  if Length(m) = 0 then
    m := Order(base);
  else
    m := m[1];
  fi;
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

fi; # end if not defined in GAP library

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
##  12274789318154414216760893584069
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# experimental helper function, uses (a+b)^p = a^p + b^p
BindGlobal("PowerModCoeffsChar",  function(coeffs, k, m, p)
  local pr, d, x, xp, xps, res, c, r, i;
  if 2*Log2Int(k) < Log2Int(p)+Length(m)+LogInt(k, p) then
    return PowerModCoeffs(coeffs, k, m);
  fi;
  pr := function(l1, l2)
    local ll;
    ll := ProductCoeffs(l1,l2);
    ReduceCoeffs(ll, m);
    ShrinkRowVector(ll);
    return ll;
  end;
  d := Length(m) - 1;
  x := [Zero(coeffs[1]), One(coeffs[1])];
  xp := PowerModCoeffs(x, p, m);
  xps := [[x[2]], xp];
  for i in [2..d-1] do
    Add(xps, pr(xps[i], xp));
  od;
  res := xps[1];
  c := ShallowCopy(coeffs);
  while k <> 0 do
    r := k mod p;
    if r = 1 then
      res := pr(res, c);
    elif r > 1 then
      res := pr(res, PowerModCoeffs(c, r, m));
    fi;
    k := (k-r)/p;
    if k <> 0 then
      # c^p
      for i in [1..Length(c)] do
        c[i] := c[i]^p;
      od;
      c := c*xps;
    fi;
  od;
  ShrinkRowVector(res);
  return res;
end);

# K finite field, cpol polynomial over K irreducible over subfield with qq
# elements and zeroes in K.
InstallGlobalFunction(FindConjugateZeroes, function(K, cpol, qq)
  local p, fac, d, r, rr, g1, g2, z, i;
  if IsPolynomial(cpol) then
    cpol := One(K)*CoefficientsOfUnivariatePolynomial(cpol);
  fi;
  if Characteristic(K) = 2 then
    return FindConjugateZeroesChar2(K, cpol, qq);
  fi;
  p := Characteristic(K);
  fac := cpol;
  d := Length(fac)-1;
  while Length(fac) > 2 do
    Info(InfoStandardFF, 4, "deg(fac) = ",Length(fac)-1, "\n");
    r := 0;
    while IsZero(r) do
      r := [Random(K),One(K)];
    od;
    rr := PowerModCoeffsChar(r, (Size(K)-1)/2, fac, p);
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
  local p, n, fac, dpm1, dlen, len, f, x, l, bas, pr, c, rr, d, aa, a, 
        nfac, g, z, i, j;
  p := Characteristic(K);
  n := DegreeOverPrimeField(K);
  # Conway polynomial as polynomial coeffs over  K
  fac := One(K)*ConwayPol(p, n);
  if p = 2 then
    return FindConjugateZeroesChar2(K, fac, 2);
  fi;
  # even divisors of p-1
  dpm1 := 2*DivisorsInt((p-1)/2);
  dlen := function()
    local i;
    i := 1;
    while i < Length(dpm1) and 3 * dpm1[i+1] < 2 * len do
      i := i+1;
    od;
    return dpm1[i];
  end;
  len := Length(fac);
  f := GF(p^n);
  # we precompute the list of x^(p^i) mod Conway polynomial,
  # i < n (using GAPs arithemetic in GF(p^n) and computing Z(p,n)^(p^i)).
  x := PrimitiveElement(f);
  l := [x];
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
    Info(InfoStandardFF, 4, "deg(fac) = ",Length(fac)-1, "\n");
    # This is a slight variant of Cantor-Zassenhaus, since the factors of 
    # fac over K are linear it is sufficient to try random polynomials 
    # of form X+c, c in K.
    # We compute (X + c) ^ ((p^n-1)/d) mod fac for some small d using
    #      (X+c)^(p^i) = X^(p^i) + c^(p^i)
    # (first summands mod fac computed above, second is computation in K)
    # and
    #      (p^n-1)/d = (1 + p + p^2 + ...+ p^(n-1)) * (p-1)/d 
    c := Random(K);
    rr := ShallowCopy(l[1]);
    rr[1] := rr[1] + c;
    for i in [2..n] do
      c := c^p;
      l[i][1] := l[i][1] + c;
      rr := pr(rr, l[i]);
      l[i][1] := l[i][1] - c;
    od;
    # we power by (p-1)/d such that rr - a (with a^d = 1) has good probability
    # to have non-trivial, and then low degree, gcd with fac
    d := dlen();
    rr := PowerModCoeffs(rr, (p-1)/d, fac);
    # now Gcd of this power - a with fac is likely to split fac 
    # (where a^d = 1)
    aa := Z(p)^((p-1)/d);
    a := aa;
    nfac := fac;
    for j in [1..Minimum(d, 5)] do
      rr[1] := rr[1] - a;
      g := GcdCoeffs(rr, fac);
      rr[1] := rr[1] + a;
      if Length(g) > 1 and Length(g) < Length(nfac) then
        nfac := g;
        if Length(nfac) = 2 then 
          break;
        fi;
      fi;
      a := a*aa;
    od;
    fac := nfac;
    if Length(fac) < len  then
      len := Length(fac);
      if len > 2 then 
        # in case of splitting we can from now compute modulo found factor fac
        for a in l do
          ReduceCoeffs(a, fac);
          ShrinkRowVector(a);
        od;
      fi;
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
    Info(InfoStandardFF, 4, "deg(fac) = ",Length(fac)-1, "\n");
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
##    "https://maths-people.anu.edu.au/~brent/ftp/factors/factors.gz",
##    false);
### newer and more factors collected by Jonathan Crombie, check occasionally
### if newer updates are available and adjust date part of URL:
##  FetchMoreFactors(
##    "http://myfactorcollection.mooo.com:8090/brentdata/May31_2022/factors.gz",
##    true);
InstallGlobalFunction(FetchMoreFactors,   function ( url, write )
  local  str, get, comm, rows, b, k, a, dir;

  # Fetch the file from R. P. Brent's ftp site and gunzip it into 'str'.

  str := "";
  get := OutputTextString(str, false);
  comm := Concatenation("wget --no-check-certificate -q ", 
                         url, " -O - | gzip -dc ");
  Process(DirectoryCurrent(), Filename(DirectoriesSystemPrograms(),"sh"),
          InputTextUser(), get, ["-c", comm]);

  rows := SplitString(str, "", "\n");
  str := 0;
  for a in rows do
    b := List(SplitString(a, "", "+- \n"), Int);
    if not IsBound(BRENTFACTORS[b[1]]) then
      BRENTFACTORS[b[1]] := [];
    elif not IsMutable(BRENTFACTORS[b[1]]) then
      BRENTFACTORS[b[1]] := ShallowCopy(BRENTFACTORS[b[1]]);
    fi;
    if '-' in a then
      k := b[2];
    else
      k := 2*b[2];
    fi;
    if not IsBound(BRENTFACTORS[b[1]][k]) then
      BRENTFACTORS[b[1]][k] := [b[3]];
    else
      if not IsMutable(BRENTFACTORS[b[1]][k]) then 
        BRENTFACTORS[b[1]][k] := ShallowCopy(BRENTFACTORS[b[1]][k]);
      fi;
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
##  <P/>
##  WARNING: The result  of this function may not be  a valid Brauer character
##  for  the  table <A>tab</A>  (that  is  an  integer linear  combination  of
##  irreducible Brauer  characters in  <A>tab</A>). For  a proper  handling of
##  several lifts  the data structure of  Brauer character tables needs  to be
##  extended (it must refer to the lift), and then the result of this function
##  should return a  Brauer character of another table that  refers to another
##  lift.
##  <Example>gap> tab := BrauerTable("M", 19);
##  BrauerTable( "M", 19 )
##  gap> # cannot translate some values to different lift
##  gap> fail in AsList(StandardValuesBrauerCharacter(tab, Irr(tab)[16]));
##  true
##  gap> # but table contains the irreducible Brauer characters for any lift
##  gap> ForAll(Irr(tab), bch-> IsGaloisInvariant(tab, bch));
##  true
##  gap> tab := BrauerTable("A18", 3);
##  BrauerTable( "A18", 3 )
##  gap> # here different lifts lead to different Brauer character tables
##  gap> bch := Irr(tab)[38];;
##  gap> IsGaloisInvariant(tab, bch);
##  false
##  gap> new := StandardValuesBrauerCharacter(tab, bch);;
##  gap> fail in AsList(new);
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


##  <#GAPDoc Label="FrobeniusCharacterValues">
##  <ManSection>
##  <Heading>Frobenius character values</Heading>
##  <Func Name="SmallestDegreeFrobeniusCharacterValue" Arg="cyc, p"/>
##  <Returns>a positive integer or <K>fail</K></Returns>
##  <Func Name="StandardFrobeniusCharacterValue" Arg="cyc, F"/>
##  <Returns>an element of <A>F</A> or <K>fail</K></Returns>
##  <Description>
##  The argument <A>cyc</A> must be a cyclotomic whose conductor and 
##  denominator are not  divisible by the prime integer <A>p</A> or 
##  the characteristic of the standard finite field <A>F</A>. 
##  <P/>
##  The order of the multiplicative group of <A>F</A> must be divisible
##  by the conductor of <A>cyc</A>. 
##  <P/>
##  Then <Ref Func="StandardFrobeniusCharacterValue"/> returns the image 
##  of <A>cyc</A> in <A>F</A> under the homomorphism which maps the
##  root of unity <C>E(n)</C> to the <Ref Func="StandardCyclicGenerator"/>
##  of order <C>n</C> in <A>F</A>. If the conditions are not fulfilled 
##  the function returns <K>fail</K>.
##  <P/>
##  The function <Ref Func="SmallestDegreeFrobeniusCharacterValue"/> returns
##  the smallest degree of a field over the prime field of order <A>p</A> 
##  containing the image of <A>cyc</A>.
##  <Example>gap> SmallestDegreeFrobeniusCharacterValue(E(13), 19);
##  12
##  gap> F := FF(19,12);
##  FF(19, 12)
##  gap> x := StandardFrobeniusCharacterValue(E(13),F);;
##  gap> x^13;
##  ZZ(19,12,[1])
##  gap> x = StandardCyclicGenerator(F, 13);
##  true
##  gap> cc := (E(13)+1/3)^4;;
##  gap> xx := StandardFrobeniusCharacterValue(cc, F);;
##  gap> xx = StandardFrobeniusCharacterValue(E(13)+1/3, F)^4;
##  true
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
# reverse direction: Frobenius character values
InstallGlobalFunction(SmallestDegreeFrobeniusCharacterValue, function(cyc, p)
  local N, c;
  N := Conductor(cyc);
  c := CoeffsCyc(cyc, N);
  if N mod p = 0 then
    return fail;
  fi;
  return OrderMod(p, N);
end);
InstallGlobalFunction(StandardFrobeniusCharacterValue, function(cyc, F)
  local p, N, c, x, res, o, j;
  p := Characteristic(F);
  N := Conductor(cyc);
  c := CoeffsCyc(cyc, N);
  if N mod p = 0 then
    return fail;
  fi;
  if (Size(F)-1) mod N <> 0 then
    return fail;
  fi;
  # image of E(N)
  x := StandardCyclicGenerator(F, N);
  res := Zero(F);
  o := Z(p)^0;
  for j in [N, N-1..1] do
    res := res * x + c[j]*o;
  od;
  return res;
end);
  
# utilities for multivariate polynomials

##  Degree of multivariate polynomial in variable (given a number or polynomial)
##  (get coefficients with PolynomialCoefficientsOfPolynomial)
InstallMethod(Degree, ["IsPolynomial", "IsObject"], 
function(pol, var)
  local ext, res, a, i, j;
  if not IsInt(var) then
    var := IndeterminateNumberOfLaurentPolynomial(var);
  fi;
  ext := ExtRepPolynomialRatFun(pol);
  res := 0;
  for i in [1,3..Length(ext)-1] do
    a := ext[i];
    for j in [1,3..Length(a)-1] do
      if a[j] = var and a[j+1] > res then
        res := a[j+1];
      fi;
    od;
  od;
  return res;
end);

# helper: numbers of indeterminates of a multivariate polynomial 
# l can also be cyclotomic or finite field element or recursive list
# of these
InstallMethod(Indets, ["IsObject"], function(l)
  local odds, is, ext, iss, res, i;
  odds := l-> l{[1,3..Length(l)-1]};
  if IsList(l) then
    is := Concatenation(List(l, Indets));
  elif IsCyc(l) or IsFFE(l) then
    is := [];
  elif IsPolynomial(l) then
    ext := ExtRepPolynomialRatFun(l);
    is := Concatenation(List(odds(ext), a-> odds(a)));
  fi;
  iss := [];
  IsSet(iss);
  res := [];
  for i in is do
    if not i in iss then
      AddSet(iss, i);
      Add(res, i);
    fi;
  od;
  return res;
end);
