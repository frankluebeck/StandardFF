#############################################################################
##
#A  TowerArith.gi                                                Frank Lübeck
##
##  The  files TowerArith.g{i,d}  contain code  to generate  and compute  in
##  towers of finite fields.
##  


##  Degree of multivariate polynomial in variable (given a number or polynomial)
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

##  <#GAPDoc Label="FiniteFieldTower">
##  <ManSection>
##  <Heading>Creating towers of finite fields</Heading>
##  <Func Name="FiniteFieldTower" Arg="K, pols..." />
##  <Func Name="ExtendFiniteFieldTower" Arg="T, pols..." />
##  <Func Name="PreviousFiniteFieldTower" Arg="T" />
##  <Func Name="BaseFiniteFieldTower" Arg="T" />
##  <Returns>a field</Returns>
##  <Description>
##  Let   <A>K</A>   be    a   field   and   <A>pols</A>   be    a   list   of
##  multivariate  polynomials  over <A>K</A>.  As  explained  in section  <Ref
##  Sect="sec:DefTower"/> it is assumed  that each polynomial contains exactly
##  one variable  that does not  occur in  previous polynomials and  that this
##  polynomial modulo  the previous  ones is  irreducible. The  first function
##  <Ref  Func="FiniteFieldTower"/>  returns  a  field which  is  a  tower  of
##  algebraic extensions,  where each  step corresponds  to one  polynomial in
##  <A>pols</A>.  It is  not checked  that the  given polynomials  fulfill the
##  required conditions. If <A>pols</A> is empty then <A>K</A> is returned.
##  <P/>
##  The second  function <Ref Func="ExtendFiniteFieldTower"/> gets  a tower of
##  extensions  <A>T</A> from  a previous  call and  extends it  by the  given
##  further polynomials <A>pols</A>.
##  <P/>
##  The  third function <Ref  Func="PreviousFiniteFieldTower"/>  returns for  a
##  tower of field  extensions the largest proper subfield in  the tower (just
##  <A>K</A> if it is a tower with only one step).
##  <P/>
##  The function <Ref Func="BaseFiniteFieldTower"/> returns the base field <A>K</A> of a field tower constructed by the previous commands.
##  <P/>
##  Despite their names these basic  constructions will probably work with any
##  available  field in  &GAP;,  but we  have  only used  and  tested it  with
##  finite  fields.  In  all  applications  in this package  we  use  it  with
##  <A>K</A> a finite prime  field <C>GF(p)</C>, see <Ref BookName="Reference"
##  Func="GF"/>.
##  
##  <P/> 
##  <Example><![CDATA[gap> x := Indeterminate(GF(13), "x");; 
##  gap> y := Indeterminate(GF(13), "y");;
##  gap> pol1 := x^2 - 2; pol2 := y^2 - x;
##  x^2+Z(13)^7
##  y^2-x
##  gap> # printed by size of base field and degrees of polynomials
##  gap> F1 := FiniteFieldTower(GF(13), pol1);
##  Tower(13;2)
##  gap> PreviousFiniteFieldTower(F1);
##  GF(13)
##  gap> F := ExtendFiniteFieldTower(F1, pol2);
##  Tower(13;2,2)
##  gap> PreviousFiniteFieldTower(F);
##  Tower(13;2)
##  gap> # elements are printed by coefficients wrt. tower basis
##  gap> gen := GeneratorsOfField(F);
##  [ <Tower(13;2,2):[ 0, 1, 0, 0 ]>, <Tower(13;2,2):[ 0, 0, 1, 0 ]> ]
##  gap> Sum(gen); Product(gen);
##  <Tower(13;2,2):[ 0, 1, 1, 0 ]>
##  <Tower(13;2,2):[ 0, 0, 0, 1 ]>
##  gap> Random(F)^(Size(F)-1);
##  <Tower(13;2,2):[ 1, 0, 0, 0 ]>
##  gap> Zero(F); One(F);
##  <Tower(13;2,2):[ 0, 0, 0, 0 ]>
##  <Tower(13;2,2):[ 1, 0, 0, 0 ]>
##  ]]></Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>

InstallGlobalFunction(FiniteFieldTower, function(F, pols...)
  local ppols, npol, prev, var, nam, coeffs, deg, res, z;
  # allow prime power or field
  if IsInt(F) then
    F := GF(F);
  fi;
  # allow list of polynomials
  if Length(pols) = 1 and IsList(pols[1]) then
    pols := pols[1];
  fi;
  # allow empty pols?
  if Length(pols) = 0 then 
    return F; 
  else
    ppols := pols{[1..Length(pols)-1]};
    npol := pols[Length(pols)];
    prev := FiniteFieldTower(F, ppols);
    var := Difference(Indets(npol), Indets(ppols))[1];
    nam := String(Indeterminate(F, var));
    coeffs := PolynomialCoefficientsOfPolynomial(npol, var);
    if IsFiniteFieldTower(prev) then
      coeffs := List(coeffs, c-> FiniteFieldTowerElement(prev, c));
    else
      # constant polynomials
      z := Zero(F);
      coeffs := List(coeffs, c-> Value(c, z));
      ConvertToVectorRep(coeffs);
    fi;
    deg := Length(coeffs)-1;
    res := rec(prev:=prev, var := var, coeffs := coeffs, nam := nam,
               deg := deg, poly := npol);
    Objectify(FiniteFieldTowerType, res);
    if HasDegreeOverPrimeField(prev) then
      SetDegreeOverPrimeField(res, DegreeOverPrimeField(prev)*deg);
    fi;
    SetCharacteristic(res, Characteristic(F));
    return res;
  fi;
end);

InstallGlobalFunction(ExtendFiniteFieldTower, function(tower, pols...)
  local npol, vars, pr, var, coeffs, nam, deg, res;
  if Length(pols) = 1 and IsList(pols[1]) then
    pols := pols[1];
  fi;
  if Length(pols) = 0 then
    return tower;
  fi;
  if not IsFiniteFieldTower(tower) then
    return FiniteFieldTower(tower, pols);
  fi;
  npol := pols[1];
  vars := [tower!.var];
  pr := tower!.prev;
  while IsFiniteFieldTower(pr) do
    Add(vars, pr!.var);
    pr := pr!.prev;
  od;
  var := Difference(Indets(npol), vars)[1];
  coeffs := PolynomialCoefficientsOfPolynomial(npol, var);
  coeffs := List(coeffs, c-> FiniteFieldTowerElement(tower, c));
  nam  := String(Indeterminate(pr, var));
  deg := Length(coeffs)-1;
  res := rec(prev:=tower, var := var, coeffs := coeffs, nam := nam,
               deg := deg, poly := npol);
  Objectify(FiniteFieldTowerType, res);
  if HasDegreeOverPrimeField(tower) then
    SetDegreeOverPrimeField(res, DegreeOverPrimeField(tower)*deg);
  fi;
  SetCharacteristic(res, Characteristic(tower));
  if Length(pols) > 1 then
    res := ExtendFiniteFieldTower(res, pols{[2..Length(pols)]});
  fi;
  return res;
end);
InstallGlobalFunction(PreviousFiniteFieldTower, function(tower)
  if IsFiniteFieldTower(tower) then
    return tower!.prev;
  else
    return Error("No previous field in finite field tower.");
  fi;
end);

InstallMethod(BaseFiniteFieldTower, ["IsFiniteFieldTower"], function(tower)
  while IsFiniteFieldTower(tower) do
    tower := tower!.prev;
  od;
  return tower;
end);

##  <#GAPDoc Label="EltPoly">
##  <ManSection>
##  <Heading>Elements by polynomials</Heading>
##  <Oper Name="FiniteFieldTowerElement" Arg="F, pol" />
##  <Returns>an element of finite field tower</Returns>
##  <Func Name="PolynomialFiniteFieldTowerElement" Arg="c" />
##  <Returns>a polynomial</Returns>
##  <Description>
##  The argument <A>F</A> must be a tower of finite fields and <A>pol</A> must
##  be  a  (multivariate)  polynomial  in  the  variables  of  the  generating
##  polynomials  of   <A>F</A>.  Then   <Ref  Oper="FiniteFieldTowerElement"/>
##  returns the residue class of that polynomial as element of <A>F</A>.
##  <P/> 
##  Vice versa, if  <A>c</A> is an element  of a finite field  tower then <Ref
##  Func="PolynomialFiniteFieldTowerElement"  />  returns the  unique  reduced
##  polynomial representing <A>c</A>, see <Ref Sect="sec:DefTower"/>.
##  <Example><![CDATA[gap> # F as above
##  gap> pol := x^4 - x^2*y^3 - 2;
##  -x^2*y^3+x^4+Z(13)^7
##  gap> c := FiniteFieldTowerElement(F, pol);
##  <Tower(13;2,2):[ 2, 0, 0, 11 ]>
##  gap> redpol := PolynomialFiniteFieldTowerElement(c);
##  Z(13)^7*x*y+Z(13)
##  gap> c = FiniteFieldTowerElement(F, redpol);
##  true]]>
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# element from multivariate polynomial
InstallMethod(FiniteFieldTowerElement, ["IsFiniteFieldTower", "IsPolynomial"],
function(fft, pol)
  local rep, z, res;
  rep := PolynomialCoefficientsOfPolynomial(pol, fft!.var);
  if IsFiniteFieldTower(fft!.prev) then
    rep := List(rep, a-> FiniteFieldTowerElement(fft!.prev, a));
  else
    # constant polys
    z := Zero(fft!.prev);
    rep := List(rep, pl-> Value(pl, z));
    if Length(rep) = 0 then
      rep := [z];
      ConvertToVectorRep(rep);
      rep := rep{[]};
    else
      ConvertToVectorRep(rep);
    fi;
  fi;
  if Length(rep) >= Length(fft!.coeffs) then
    if not IsMutable(rep) then 
      rep := ShallowCopy(rep);
    fi;
    ReduceCoeffs(rep, fft!.coeffs);
    ShrinkRowVector(rep);
  fi;
  res := rec(field := fft, rep := rep);
  Objectify(FiniteFieldTowerElementType, res);
  return res;
end);
# prime field with constant polynomial
InstallOtherMethod(FiniteFieldTowerElement, ["IsFinite and IsFFECollection", "IsPolynomial"],
function(gf, pl)
  if Length(Indets(pl)) > 0 then
    Error("polynomial must be constant");
  fi;
  return Value(pl, Zero(gf));
end);

# element from coeff list wrt. previous field
# ... to be documented for the user ???
InstallMethod(FiniteFieldTowerElement, ["IsFiniteFieldTower", "IsList"], 
function(K, c)
  local res;
  c := ShallowCopy(c);
  while Length(c)>0 and  IsZero(c[Length(c)]) do
    Remove(c);
  od;
  if not IsFiniteFieldTower(K!.prev) then
    ConvertToVectorRep(c);
  fi;
  if Length(c) > K!.deg then
    ReduceCoeffs(c, K!.coeffs);
  fi;
  res := rec(field := K, rep := c);
  Objectify(FiniteFieldTowerElementType, res);
  return res;
end);

InstallGlobalFunction(PolynomialFiniteFieldTowerElement, function(a)
  local L, base, v, rep, len, res, i;
  if IsFFE(a) then
    return a * Indeterminate(GF(Characteristic(a)))^0;
  fi;
  L := a!.field;
  base := BaseFiniteFieldTower(L);
  v := Indeterminate(base, L!.var);
  rep := a!.rep;
  len := Length(rep);
  if len = 0 then
    return 0*v;
  else
    if IsFiniteFieldTowerElement(rep[1]) then
      rep := List(rep, PolynomialFiniteFieldTowerElement);
    fi;
    res := rep[len]*v^0;
    for i in [len-1,len-2..1] do
      res := res*v+rep[i];
    od;
    return res;
  fi;
end);

##  <#GAPDoc Label="EltVec">
##  <ManSection>
##  <Heading>Elements with respect to tower basis</Heading>
##  <Meth Name="AsVector" Label="for elements in field towers" Arg="c"/>
##  <Returns>a vector over the base field</Returns>
##  <Func Name="FiniteFieldTowerElementVector" Arg="F, vec"/>
##  <Returns>an element of a finite field tower</Returns>
##  <Description>
##  If  <A>c</A> is  an element  of a  finite field  tower <C>F</C>  then <Ref
##  Meth="AsVector"  Label="for  elements  in  field towers"  />  returns  the
##  coordinate vector with respect to the  tower basis of <C>F</C> (entries in
##  the base field of the tower), see <Ref Sect="sec:DefTower"/>.
##  <P/> 
##  The reverse map is provided by <Ref Func="FiniteFieldTowerElementVector"/>.
##  <Example><![CDATA[gap> # F as above
##  gap> c := Sum(GeneratorsOfField(F))^123456;
##  <Tower(13;2,2):[ 6, 7, 0, 10 ]>
##  gap> vec := AsVector(c);
##  [ Z(13)^5, Z(13)^11, 0*Z(13), Z(13)^10 ]
##  gap> c = FiniteFieldTowerElementVector(F, vec);
##  true]]>
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# coefficients over base field wrt. tower basis
InstallMethod(AsVector, ["IsFiniteFieldTowerElement"], function(x)
  local L, rep, i;
  L := x!.field;
  rep := x!.rep;
  if Length(rep) < L!.deg then
    rep := ShallowCopy(rep);
    for i in [Length(rep)+1..L!.deg] do
      Add(rep, Zero(L!.prev));
    od;
  fi;
  if IsFiniteFieldTower(L!.prev) then
    return Concatenation(List(rep, AsVector));
  else
    return rep;
  fi;
end);
# inverse
InstallGlobalFunction(FiniteFieldTowerElementVector,  function(K, v)
  local prev, m, rep, res;
  prev := K!.prev;
  if IsFiniteFieldTower(prev) then
    m := Length(v)/K!.deg;
    rep := List([0..K!.deg-1], i-> FiniteFieldTowerElementVector(prev,
              v{[i*m+1..(i+1)*m]}));
  else
    rep := ShallowCopy(v);
  fi;
  while Length(rep)>0 and IsZero(rep[Length(rep)]) do
    Unbind(rep[Length(rep)]);
  od;
  res := rec(field := K, rep := rep);
  Objectify(FiniteFieldTowerElementType, res);
  return res;
end);

##  <#GAPDoc Label="EltSteinitz">
##  <ManSection>
##  <Heading>Elements by Steinitz number</Heading>
##  <Oper Name="SteinitzNumber" Arg="c"/>
##  <Returns>an integer</Returns>
##  <Oper Name="ElementSteinitzNumber" Arg="F, nr"/>
##  <Returns>an element in finite field tower</Returns>
##  <Description>
##  Here  <A>c</A> must  be  an  element of  a  finite  field tower  <A>F</A>.
##  These  functions only  work if  the  base field  of <A>F</A>  is a  finite
##  prime  field   <C>GF(p)</C>  (or   <C>FF(p,1)</C>),  see   <Ref  Func="GF"
##  BookName="Reference"/>.
##  <P/> 
##  The first function returns the Steinitz  number of the element <A>c</A> as
##  defined in&nbsp;<Ref Sect="sec:DefTower"/>. And the second function is the
##  reverse and returns the field element with given Steinitz number.
##  <Example><![CDATA[gap> # F as above, base field is GF(13)
##  gap> List(GeneratorsOfField(F), SteinitzNumber);
##  [ 13, 169 ]
##  gap> ElementSteinitzNumber(F, Size(F)-1);
##  <Tower(13;2,2):[ 12, 12, 12, 12 ]>]]>
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
# only works if base field of tower is prime field GF(p)
InstallMethod(SteinitzNumber, ["IsFFE"], IntFFE);
InstallMethod(SteinitzNumber, ["IsFiniteFieldTowerElement"], function(x)
  local F, q, rep;
  F := x!.field;
  q := Size(F!.prev);
  rep := x!.rep;
  if Length(rep) = 0 then
    return 0;
  else
    return ValuePol(List(rep, SteinitzNumber), q);
  fi;
end);
# element from Steinitz number
InstallOtherMethod(ElementSteinitzNumber, 
                               ["IsStandardPrimeField", "IsInt"],
function(K, nr)
  return nr*One(K);
end);
InstallMethod(ElementSteinitzNumber, ["IsFiniteFieldTower", "IsInt"],
function(K, nr)
  local v, p, d;
  p := Characteristic(K);
  d := DegreeOverPrimeField(K);
  v := CoefficientsQadic(nr, p);
  while Length(v) < d do
    Add(v, 0);
  od;
  v := One(GF(p))*v;
  ConvertToVectorRep(v, p);
  return FiniteFieldTowerElementVector(K, v);
end);


InstallOtherMethod(ZeroOp, ["IsFiniteFieldTower"], function(fft)
  local c;
  if not IsFiniteFieldTower(fft!.prev) then
    c := [Zero(fft!.prev)];
    ConvertToVectorRep(c);
    c := c{[]};
  else
    c := [];
  fi;
  return FiniteFieldTowerElement(fft, c);
end);
InstallMethod(Zero, ["IsFiniteFieldTower"], function(fft)
  return ZeroMutable(fft);
end);
InstallOtherMethod(OneOp, ["IsFiniteFieldTower"], function(fft)
  local c;
  c := [One(fft!.prev)];
  if IsFFE(c[1]) then
    ConvertToVectorRep(c);
  fi;
  return FiniteFieldTowerElement(fft, c);
end);
InstallMethod(One, ["IsFiniteFieldTower"], fft-> OneOp(fft));

BindGlobal("LastGenerator", function(fft)
  local c;
  if not IsFiniteFieldTower(fft) then
    Error("no field tower");
  else
    c := [Zero(fft!.prev), One(fft!.prev)];
    return FiniteFieldTowerElement(fft, c);
  fi;
end);

InstallMethod(GeneratorsOfField, ["IsFiniteFieldTower"], function(fft)
  local pgen;
  if IsFiniteFieldTower(fft!.prev) or (not IsPrimeField(fft!.prev)) then
    pgen := ShallowCopy(GeneratorsOfField(fft!.prev));
    pgen := List(pgen, a-> FiniteFieldTowerElement(fft, [a]));
    Add(pgen, LastGenerator(fft));
    return pgen;
  else
    return [LastGenerator(fft)];
  fi;
end);

InstallMethod(ViewString, ["IsFiniteFieldTower"], function(fft)
  local res;
  if not IsFiniteFieldTower(fft!.prev) then
    return Concatenation("Tower(", String(Size(fft!.prev)),
                                ";",String(fft!.deg),")");
  fi;
  res := ShallowCopy(ViewString(fft!.prev));
  Remove(res);
  Append(res,",");
  Append(res, String(fft!.deg));
  Append(res,")");
  return res;
end);
InstallMethod(ViewObj, ["IsFiniteFieldTower"], 
                       function(fft) Print(ViewString(fft)); end);
InstallMethod(ViewString, ["Is8BitVectorRep"], function(v)
  if Length(v) = 0 then
    return "[]";
  else
    TryNextMethod();
  fi;
end);
InstallMethod(ViewString, ["IsFiniteFieldTowerElement"], function(x)
  local l;
  l := AsPlist(AsVector(x));
  if IsPrimeField(BaseFiniteFieldTower(x!.field)) then
    l := List(l, IntFFE);
  fi;
  return Concatenation("<", ViewString(x!.field), ":", 
                        ViewString(l), ">");
end);

InstallMethod(Size, ["IsFiniteFieldTower"], function(fft) 
  return Size(fft!.prev^fft!.deg);
end);


InstallMethod(IsZero, ["IsFiniteFieldTowerElement"], function(x)
  return ForAll(x!.rep, IsZero);
end);
InstallMethod(ZeroOp, ["IsFiniteFieldTowerElement"], function(x)
  return Zero(x!.field);
end);
InstallMethod(IsOne, ["IsFiniteFieldTowerElement"], function(x)
  return Length(x!.rep) > 0 and IsOne(x!.rep[1]) and
                            ForAll([2..Length(x!.rep)], i-> IsZero(x!.rep[i]));
end);
InstallMethod(OneOp, ["IsFiniteFieldTowerElement"], function(x)
  return One(x!.field);
end);

# arithetic and random is easily implemented by recursion
InstallMethod(\+, ["IsFiniteFieldTowerElement", "IsFiniteFieldTowerElement"], 
function(x, y)
  if IsIdenticalObj(x!.field, y!.field) then
    return FiniteFieldTowerElement(x!.field, x!.rep + y!.rep);
  else
    return fail;
  fi;
end);

# helper for multiplication, we cache reductions of var^i for
# deg <= i < 2*deg
InstallMethod(ReducedGeneratorPowers, ["IsFiniteFieldTower"], 
function(fft)
  local z, zl, o, deg, res, v, vv, i;
  z := Zero(fft!.prev);
  zl := [z];
  if IsFFE(z) then
    ConvertToVectorRep(zl);
  fi;
  o := One(fft!.prev);
  deg := fft!.deg;
  res := [];
  v := List([1..deg-1], i-> z);
  Add(v, o);
  if IsFFE(z) then
    ConvertToVectorRep(v);
  fi;
  for i in [deg+1..2*deg-1] do
    vv := Concatenation(zl,v);
    ReduceCoeffs(vv, fft!.coeffs);
    while Length(vv) < deg do
      Append(vv, zl);
    od;
    Add(res, vv{[1..deg]});
    v := res[i-deg];
  od;
  return res;
end);

InstallMethod(\*, ["IsFiniteFieldTowerElement", "IsFiniteFieldTowerElement"], 
function(x, y)
  local K, deg, res;
  K := x!.field;
  deg := K!.deg;
  if Length(x!.rep) = 0 then
    return x;
  elif Length(y!.rep) = 0 then
    return y;
  fi;
  if IsIdenticalObj(K, y!.field) then
    res := ProductCoeffs(x!.rep, y!.rep);
    if Length(res) > deg then
      res := res{[1..deg]} 
                       + res{[deg+1..Length(res)]}*ReducedGeneratorPowers(K);
    fi;
    return FiniteFieldTowerElement(x!.field, res);
  else
    return fail;
  fi;
end);
##  # sparse version
##  InstallMethod(\*, ["IsFiniteFieldTowerElement", 
##                     "IsFiniteFieldTowerElement"], function(x, y)
##    local K, deg, si, res, k, c, ii, a;
##    K := x!.field;
##    deg := K!.deg;
##    if Length(x!.rep) = 0 then
##      return x;
##    elif Length(y!.rep) = 0 then
##      return y;
##    fi;
##    # [e, list of [i, -a_i]]
##    si := K!.sparseinfo;
##    if IsIdenticalObj(K, y!.field) then
##      res := ProductCoeffs(x!.rep, y!.rep);
##      k := Length(res);
##      while k > deg do
##        c := res[k];
##        if not IsZero(c) then
##          for a in si do
##            ii := k-deg+a[1];
##            res[ii] := res[ii]+c*a[2];
##          od;
##        fi;
##        k := k-1;
##      od;
##      if Length(res) > deg then
##        res := res{[1..deg]};
##      fi;
##      return FiniteFieldTowerElement(x!.field, res);
##    else
##      return fail;
##    fi;
##  end);

InstallMethod(\*, ["IsFiniteFieldTowerElement", "IsFFE"], function(x, c)
  return FiniteFieldTowerElement(x!.field, c * x!.rep);
end);
InstallMethod(\*, ["IsFFE", "IsFiniteFieldTowerElement"], function(c, x)
  return FiniteFieldTowerElement(x!.field, c * x!.rep);
end);
InstallMethod(\*, ["IsFiniteFieldTowerElement", "IsInt"], function(x, c)
  return FiniteFieldTowerElement(x!.field, (c*One(x!.field!.prev)) * x!.rep);
end);
InstallMethod(\*, ["IsInt", "IsFiniteFieldTowerElement"], function(c, x)
  return FiniteFieldTowerElement(x!.field, (c*One(x!.field!.prev)) * x!.rep);
end);
InstallMethod(AdditiveInverseOp, ["IsFiniteFieldTowerElement"], function(x)
  return FiniteFieldTowerElement(x!.field, -x!.rep);
end);
##  # very easy to implement, but not very efficient:
##  InstallMethod(InverseOp, ["IsFiniteFieldTowerElement"], function(x)
##    if IsZero(x) then return fail; fi;
##    return x^(Size(x!.field)-2);
##  end);
InstallMethod(InverseOp, ["IsFiniteFieldTowerElement"], function(x)
  local L, rm, r, tm, t, qr, q, tp, tt;
  if IsZero(x) then return fail; fi;
  L := x!.field;
  rm := L!.coeffs;
  r := x!.rep;
  tm := Zero(L);
  t := One(L);
  while Length(r) > 1 do
    qr := QuotRemPolList(rm, r);
    q := FiniteFieldTowerElement(L, qr[1]);
    tp := tm - q*t;
    rm := r;
    r := qr[2];
    while IsZero(r[Length(r)]) do
      Remove(r);
    od;
    tm := t;
    t := tp;
  od;
  tt := [INV(r[1])];
  if IsFFE(tt[1]) then
    ConvertToVectorRep(tt);
  fi;
  return t*FiniteFieldTowerElement(L, tt);
end);


InstallMethod(\=, ["IsFiniteFieldTowerElement", "IsFiniteFieldTowerElement"], 
function(x, y)
  return x!.field = y!.field and x!.rep = y!.rep;
end);

InstallMethod(\in, ["IsFiniteFieldTowerElement", "IsFiniteFieldTower"], 
function(x, K)
  return IsIdenticalObj(x!.field, K);
end);

InstallMethod(Random, ["IsFiniteFieldTower"], function(K)
  local c;
  c := List([1..K!.deg], i-> Random(K!.prev));
  while Length(c)>0 and IsZero(c[Length(c)]) do
    Remove(c);
  od;
  if not IsFiniteFieldTower(K!.prev) then
    ConvertToVectorRep(c);
  fi;
  return FiniteFieldTowerElement(K, c);
end);

##  <#GAPDoc Label="TowerBasis">
##  <ManSection>
##  <Attr Name="TowerBasis" Arg="F" />
##  <Returns>list of elements</Returns>
##  <Description>
##  The  argument <A>F</A>  must  be  a tower  of  finite  fields as  returned
##  by  <Ref  Oper="FiniteFieldTower"/>.  The  return value  is  the  list  of
##  elements in <A>F</A>  which form the tower basis of  <A>F</A> described in
##  section&nbsp;<Ref Sect="sec:DefTower"/>.
##  <P/> 
##  <Example><![CDATA[gap> # with F as above
##  gap> TowerBasis(F);
##  [ <Tower(13;2,2):[ 1, 0, 0, 0 ]>, <Tower(13;2,2):[ 0, 1, 0, 0 ]>, 
##    <Tower(13;2,2):[ 0, 0, 1, 0 ]>, <Tower(13;2,2):[ 0, 0, 0, 1 ]> ]]]>
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallMethod(TowerBasis, ["IsFiniteFieldTower"], function(K)
  local pbas, z, bas, rep, i, b;
  if IsFiniteFieldTower(K!.prev) then
    pbas := TowerBasis(K!.prev);
  else
    pbas := [Z(Characteristic(K))^0];
  fi;
  z := Zero(K!.prev);
  bas := [];
  for i in [1..K!.deg] do
    rep := List([1..i], j-> z);
    for b in pbas do
      rep[i] := b;
      Add(bas, FiniteFieldTowerElement(K, ShallowCopy(rep)));
    od;
  od;
  return bas;
end);

# minimal polynomial over base field
InstallOtherMethod(MinimalPolynomial, 
                           ["IsFiniteFieldTowerElement", "IsPosInt"], 
function(a, nr)
  local L, base, p, bas, mat, aa, st;
  L := a!.field;
  base := BaseFiniteFieldTower(L);
  p := Characteristic(L);
  
  # over GF(p) we use matrix wrt tower basis
  if base = GF(p) then
    bas := TowerBasis(L);
    mat := List(bas, b-> AsVector(b*a));
    ConvertToMatrixRep(mat, p);
    return MinimalPolynomial(GF(p), mat, nr);
  fi;

  # otherwise we compute powers and echelonize
  aa := One(L);
  st := FindLinearCombination(AsVector(aa), false);
  while st[2] = fail do
    aa := aa*a;
    st := FindLinearCombination(AsVector(aa), st[1]);
  od;
  return UnivariatePolynomial(base, Concatenation(-st[2], [One(base)]), nr);
end);

InstallMethod(Enumerator, ["IsFiniteFieldTower"], function(K)
  if not IsStandardPrimeField(BaseFiniteFieldTower(K)) then
    Error("Enumerator only available if base field is a GF(p).");
  fi;
  return EnumeratorByFunctions(K, rec(field := K,
    NumberElement := function(enum, x)
      return SteinitzNumber(x)+1;
    end,
    ElementNumber := function(enum, nr)
      return ElementSteinitzNumber(enum!.field, nr-1);
    end,
    Length := function(enum)
      return Size(enum!.field);
    end,
    AsList := function(enum)
      return List([0..Size(K)-1], i-> ElementSteinitzNumber(K,i));
    end));
end);
InstallMethod(Iterator, ["IsFiniteFieldTower"], 
              T-> IteratorList(Enumerator(T)));
InstallMethod(AsList, ["IsFiniteFieldTower"], 
              T-> AsList(Enumerator(T)));


