#############################################################################
##
#A  StandardFF.gi                                                Frank Lübeck
##  
##  The files StandardFF.g{i,d}  contain code to compute  standard models of
##  finite fields.  Fields are  returned as field  towers with  prime degree
##  indices, and as single extensions over the prime field.
##  

# args: p, r, k
# Let f = X^r + g(X) our standard irreducible polynomial for the 
# degree r extension over the field FF(p,r^(k-1)). This function
# returns the Steinitz number of the polynomial g(X).
InstallGlobalFunction(SteinitzNumberForPrimeDegree, function(p, r, k)
  local Fp, stpd, stpdr, q, F, o, i, nr, x, a, l;
  Fp := StandardFiniteField(p, 1);
  if not IsBound(Fp!.SteinitzPrimeDeg) then
    Fp!.SteinitzPrimeDeg := [];
  fi;
  stpd := Fp!.SteinitzPrimeDeg;
  if not IsBound(stpd[r]) then
    stpd[r] := [];
  fi;
  stpdr := stpd[r];
  if not IsBound(stpdr[k]) then
    if p = r then
      # Artin-Schreier case
      # Xk^p - Xk - \prod_{1..k-1} x_i^(p-1)
      # for q = p^(p^(k-1)) this yields (p-1)*(q + q/p)
      q := p^(p^(k-1));
      stpdr[k] := (p-1)*(q+q/p);
    elif r = 2 and p mod 4 = 3 then
      # special case for k = 1, 2
      if k = 1 then
        # X^2 + 1
        stpdr[k] := 1;
      elif k = 2 then
        # find a non-square
        F := StandardFiniteField(p, 2);
        o := One(F);
        i := 1;
        nr := StandardAffineShift(p^2, i);
        x := ElementSteinitzNumber(F, nr);
        while nr = 0 or x^((p^2-1)/2) = o do
          i := i+1;
          nr := StandardAffineShift(p^2, i);
          x := ElementSteinitzNumber(F, nr);
        od;
        # X^2 - x
        stpdr[k] := SteinitzNumber(-x);
      else
        # Xk - x{k-1}
        stpdr[k] := (p-1)*p^(2^(k-2));
      fi;
    elif (p-1) mod r = 0 then
      if k = 1 then
        # find a non r-th power
        i := 1;
        nr := StandardAffineShift(p, i);
        while nr = 0 or PowerMod(nr, (p-1)/r, p) = 1 do
          i := i+1;
          nr := StandardAffineShift(p, i);
        od;
        # X^r - nr
        stpdr[k] := p - nr;
      else
        # Xk^r - x{k-1}
        stpdr[k] := (p-1)*p^(r^(k-2));
      fi;
    else
      # in general we use pseudo random polynomials
      F := StandardFiniteField(p, r^(k-1));
      if k = 1 then
        # new generator will have norm 1
        a := -One(F);
      else
        # new generator will have previous one as norm
        a := -PrimitiveElement(F);
      fi;
      l := StandardIrreducibleCoeffList(F, r, a);
      Remove(l);
      while IsZero(l[Length(l)]) do
        Remove(l);
      od;
      stpdr[k] := ValuePol(List(l, SteinitzNumber), Size(F));
    fi;
  fi;
  return stpdr[k];
end);
# a version that tries to call external program using NTL
InstallGlobalFunction(SteinitzNumberForPrimeDegreeNTL, function(p, r, k)
  local Fp, stpd, stpdr, d, prg, inp, res;

  # for cache same as above
  Fp := StandardFiniteField(p, 1);
  if not IsBound(Fp!.SteinitzPrimeDeg) then
    Fp!.SteinitzPrimeDeg := [];
  fi;
  stpd := Fp!.SteinitzPrimeDeg;
  if not IsBound(stpd[r]) then
    stpd[r] := [];
  fi;
  stpdr := stpd[r];
  if not IsBound(stpdr[k]) then
    d := DirectoriesPackageLibrary("StandardFF","ntl");
    # currently we only have standalone programs for k = 1 and p < 2^50
    if k = 1 then
      if p = 2 then
        if r <> p and (p-1) mod r <> 0 then
          prg := Filename(d, "findstdirrGF2");
          if prg <> fail then
            inp := String(r);
            res := IO_PipeThrough(prg,[],inp);
            stpdr[k] := Int(Filtered(res, IsDigitChar));
          fi;
        fi;
      elif p < 2^50 then
        if r <> p and (p-1) mod r <> 0 then
          prg := Filename(d, "findstdirrGFp");
          if prg <> fail then
            inp := Concatenation(String(p), "\n", String(r));
            res := IO_PipeThrough(prg,[],inp);
            stpdr[k] := Int(Filtered(res, IsDigitChar));
          fi;
        fi;
      fi;
    fi;
    if not IsBound(stpdr[k]) then
      # use the GAP function
      SteinitzNumberForPrimeDegree(p, r, k);
    fi;
  fi;
  return stpdr[k];
end);

# for displaying the corresponding polynomial
InstallGlobalFunction(StandardPrimeDegreePolynomial, function(p, r, k)
  local st, qq, K, c, v;
  st := SteinitzNumberForPrimeDegree(p, r, k);
  qq := p^(r^(k-1));
  K := FF(p, r^(k-1));
  c := CoefficientsQadic(st, qq);
  while Length(c) < r do
    Add(c, 0);
  od;
  Add(c, 1);
  c := List(c, i-> AsPolynomial(ElementSteinitzNumber(K, i)));
  v := Indeterminate(FF(p,1), Concatenation("x",String(r),"_",String(k)));
  return ValuePol(c, v);
end);


# args: K, deg, lcoeffs, b
# Let f = poly(lcoeffs) + X^deg be irreducible in K[X].
# Return K[X] / f with primitve element bX mod f, where b in K
# such that bX mod f is a generator over the prime field.
# We assume that K = F[Y] / g and K knows powers of the primitive element
# Y mod g in its tower basis.  
# The returned fields also stores the powers of bX mod f in its tower
# basis.
# (The semiechelon code is similar to 'Matrix_OrderPolynomialInner'.)
InstallGlobalFunction(_ExtensionWithTowerBasis, function(K, deg, lcoeffs, b)
  local dK, zK, zKl, co, d, F, zero, one, pK, vec, vecs, pols, 
        zeroes, pmat, v, w, p, piv, x, c, vnam, var, ivar, Kp, fam, i, j;
  dK := DegreeOverPrimeField(K);
  d := dK * deg;
  F := LeftActingDomain(K);
  # shortcut for prime degree over prime field
  if dK = 1 then
    p := ShallowCopy(lcoeffs);
    while Length(p) < deg do
      Add(p, 0*p[1]);
    od;
    Add(p, One(p[1]));
    MakeImmutable(p);
    pmat := IdentityMat(deg, F);
  else
    zK := Zero(K);
    zero := Zero(F);
    one := One(F);
    pK := PrimitivePowersInTowerBasis(K);
    co := function(x)
      local res;
      if dK = 1 then
        res := [x];
      else
        res := x![1];
        if not IsList(res) then
          res := zero * [1..dK];
          res[1] := x![1];
        fi;
      fi;
      ConvertToVectorRep(res);
      return res;
    end;

    # We collect (bX)^i mod f, i = 0..d*dK-1, and compute
    # the minimal polynomial of bX over F.
    # one
    vec := NullMat(1, d, F)[1];
    vec[1] := one;
    vecs := [];
    pols := [];
    zeroes := [];
    pmat := [];
    # X
    v := NullMat(1,deg,K)[1];
    v[1] := One(K);
    for i in [1..d+1] do
      if i <= d then
        Add(pmat,vec);
      fi;
      w := ShallowCopy(vec);
      p := ShallowCopy(zeroes);
      Add(p, one);
      ConvertToVectorRepNC(p, F);
      piv := PositionNonZero(w, 0);
      # reduce
      while piv <= d and IsBound(vecs[piv]) do
          x := -w[piv];
          if IsBound(pols[piv]) then
              AddCoeffs(p, pols[piv], x);
          fi;
          AddRowVector(w, vecs[piv],  x, piv, d);
          piv := PositionNonZero(w,piv);
      od;
      if i <= d then
        x := Inverse(w[piv]);
        MultVector(p, x);
        MakeImmutable(p);
        pols[piv] := p;
        MultVector(w, x );
        MakeImmutable(w);
        vecs[piv] := w;
        Add(zeroes,zero);

        # multiply by  bX and find next vec in tower basis
        v := b*v;
        c := -v[deg];
        for j in [deg, deg-1..2] do
          v[j] := v[j-1];
        od;
        v[1] := zK;
        if not IsZero(c) then
          for j in [1..Length(lcoeffs)] do
            v[j] := v[j] + c*lcoeffs[j];
          od;
        fi;
        vec := [];
        for c in v do
          Append(vec, co(c) * pK);
        od;
        ConvertToVectorRepNC(vec, F);
      fi;
    od;
    # Now the last p is the minimal polynomial over F
    # and pmat are the primitive powers in the tower basis for the new
    # extension.
    MakeImmutable(p);
    ConvertToMatrixRep(pmat);
  fi;
  # generate the new extension
  vnam := Concatenation("x", String(d));
  var := Indeterminate(F, vnam);
  ivar := IndeterminateNumberOfUnivariateLaurentPolynomial(var);
  Kp := AlgebraicExtensionNC(F, UnivariatePolynomial(F, p, ivar), vnam);

# temporary work around for a bug in GAP-dev
if Characteristic(K) > 2^16 then
  fam := ElementsFamily(FamilyObj(Kp));;
  fam!.reductionMat := List(fam!.reductionMat,List);
fi;

  Setter(IsStandardFiniteField)(Kp, true);
  Setter(PrimitivePowersInTowerBasis)(Kp, pmat);
  if dK = 1 then
    # identity matrix
    Setter(TowerBasis)(Kp, pmat);
  fi;
  # let elements know to be in standard field
  fam := FamilyObj(RootOfDefiningPolynomial(Kp));
  fam!.extType := Subtype(fam!.extType, IsStandardFiniteFieldElement);
  fam!.baseType := Subtype(fam!.baseType, IsStandardFiniteFieldElement);
  SetFilterObj(OneImmutable(Kp), IsStandardFiniteFieldElement);
  SetFilterObj(ZeroImmutable(Kp), IsStandardFiniteFieldElement);
  SetFilterObj(RootOfDefiningPolynomial(Kp), IsStandardFiniteFieldElement);

  return Kp;
end);

##  <#GAPDoc Label="FF">
##  <ManSection>
##  <Heading>Constructing standard finite fields</Heading>
##  <Func Name="StandardFiniteField" Arg="p, n"/>
##  <Func Name="FF" Arg="p, n"/>
##  <Returns>a finite field</Returns>
##  <Func Name="StandardPrimeDegreePolynomial" Arg="p, r, k" />
##  <Returns>a polynomial of degree <A>r</A></Returns>
##  <Description>
##  The   arguments   are   a   prime  <A>p</A>   and   a   positive   integer
##  <A>n</A>.   The   function  <Ref   Func="FF"/>   (or   its  synomym   <Ref
##  Func="StandardFiniteField"/>)  is  one  of  the  main  functions  of  this
##  package.   It   returns   a   standardized   field   <C>F</C>   of   order
##  <M><A>p</A>^{<A>n</A>}</M>. It  is  implemented  as  a   simple  extension
##  over  the   prime  field  <C>GF(p)</C>  using   <Ref  BookName="Reference"
##  Oper="AlgebraicExtension" />
##  <P/>
##  The   polynomials    used   for   the   prime    degree   extensions   are
##  accessible    with   <Ref    Func="StandardPrimeDegreePolynomial"/>.   For
##  arguments  <A>p,  r,  k</A>  it  returns  the  irreducible  polynomial  of
##  degree  <A>r</A>   for  the  <A>k</A>-th  iterated   extension  of  degree
##  <A>r</A>  over  the  prime  field.  The  polynomial  is  in  the  variable
##  <C>x</C><Emph>r</Emph><C>_</C><Emph>k</Emph>  and   the  coefficients  can
##  contain  variables <C>x</C><Emph>r</Emph><C>_</C><Emph>l</Emph>  with <M>l
##  &lt; k</M>.
##  <P/>
##  <Example>gap> Fp := FF(2, 1);
##  GF(2)
##  gap> F := FF(2, 100);
##  FF(2, 100)
##  gap> Size(F);
##  1267650600228229401496703205376
##  gap> p := NextPrimeInt(10^50);
##  100000000000000000000000000000000000000000000000151
##  gap> K := FF(p, 60);
##  FF(100000000000000000000000000000000000000000000000151, 60)
##  gap> LogInt(Size(K), 10);
##  3000
##  gap> F := FF(13, 9*5);
##  FF(13, 45)
##  gap> StandardPrimeDegreePolynomial(13, 3, 1);
##  x3_1^3+Z(13)^7
##  gap> StandardPrimeDegreePolynomial(13, 3, 2);
##  x3_2^3-x3_1
##  gap> StandardPrimeDegreePolynomial(13, 5, 1);
##  x5_1^5+Z(13)^4*x5_1^2+Z(13)^4*x5_1-Z(13)^0
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallGlobalFunction(StandardFiniteField, function(p, n)
  local F, id, ext, fac, lf, n1, nK, st, q1, l, K, lK, c, L, b;
  if not IsPrimeInt(p) then
    Error("StandardFiniteField: first argument must be a prime\n");
  fi;
  F := GF(p);
  # we cache extensions in prime field
  if not IsBound(F!.extensions) then
    SetIsStandardPrimeField(F, true);
    F!.extensions := [F];
    id := IdentityMat(1, F);
    ConvertToMatrixRep(id, p);
    Setter(PrimitivePowersInTowerBasis)(F, id);
  fi;
  ext := F!.extensions;
  if IsBound(ext[n]) then
    return ext[n];
  fi;
  
  # construct field recursively via prime degree extensions
  fac := Collected(Factors(n));
  lf := fac[Length(fac)];
  n1 := lf[1]^(lf[2]-1);
  nK := n/lf[1];
  st := SteinitzNumberForPrimeDegree(p, lf[1], lf[2]);
  q1 := p^n1;
  l := CoefficientsQadic(st, q1);
  K := StandardFiniteField(p, nK);
  lK := List(l, y-> EmbedSteinitz(p, n1, nK, y));
  c := List(lK, nr-> ElementSteinitzNumber(K, nr));
  # primitive element of new extension is product of prime power
  # degree generators, class of bX, we construct b (product of the 
  # other prime power degree generators):
  b := ElementSteinitzNumber(K, p^(Position(StdMonDegs(nK), nK/n1)-1));
  L := _ExtensionWithTowerBasis(K, lf[1], c, b);
  ext[n] := L;
  return L;
end);

# default for other fields
InstallOtherMethod(IsStandardFiniteField, ["IsField"], ReturnFalse);

# maps between elements of standard finite field and its
# standard finite field tower
InstallMethod(TowerBasis, ["IsStandardFiniteField"], function(F)
  return PrimitivePowersInTowerBasis(F)^-1;
end);

# elements as coefficient vectors in tower basis
InstallMethod(AsVector, ["IsStandardFiniteFieldElement"], function(x)
  local F;
  F := FamilyObj(x)!.wholeField;
  return ExtRepOfObj(x) * PrimitivePowersInTowerBasis(F);
end);

InstallOtherMethod(AsVector, 
    ["IsStandardFiniteField", "IsStandardFiniteFieldElement"], 
function(F, x)
  return ExtRepOfObj(x) * PrimitivePowersInTowerBasis(F);
end);

# elements as polynomials
InstallMethod(GeneratorMonomials, [IsStandardFiniteField], 
function(F)
  local res, d, f, Fp, s, a, i, v;
  res := rec(vars := []);
  d := DegreeOverPrimeField(F);
  if d = 1 then
    return res;
  fi;
  f := Collected(Factors(d));
  Fp := PrimeField(F);
  for a in f do
    res.(a[1]) := rec();
    for i in [1..a[2]] do
      s := Concatenation("x",String(a[1]),"_",String(i));
      v := Indeterminate(Fp, s);
      res.(a[1]).(i) := v;
      Add(res.vars, v);
    od;
  od;
  return res;
end);
InstallMethod(TowerBasisMonomials, [IsStandardFiniteField],
function(F)
  local d, m, l, o, res, b, i, a;
  d := DegreeOverPrimeField(F);
  m := GeneratorMonomials(F);
  l := StdMon(d)[1];
  o := One(PrimeField(F));
  res := [];
  for a in l do
    b := o;
    i := 1;
    while i < Length(a) do
      b := b*m.(a[i]).(a[i+1])^a[i+2];
      i := i+3;
    od;
    Add(res, b);
  od;
  return res;
end);
InstallOtherMethod(AsPolynomial, ["IsFFE"], IdFunc);
InstallMethod(AsPolynomial, ["IsStandardFiniteFieldElement"],
function(x)
  local F;
  F := FamilyObj(x)!.wholeField;
  return AsVector(x) * TowerBasisMonomials(F);
end);
# vice versa, we also support non-reduced polynomials
InstallMethod(ElementPolynomial, [IsStandardFiniteField, IsPolynomial],
function(F, pol)
  local gm, tbm, poss, a, tb, res;
  gm := GeneratorMonomials(F);
  if not IsBound(F!.gmeltvars) then
    tbm := TowerBasisMonomials(F);
    poss := List(gm.vars, x-> Position(tbm, x));
    a := PrimitiveElement(F);
    tb := TowerBasis(F);
    F!.gmeltvars := List(poss, i-> ValuePol(tb[i], a)); 
  fi;
  res := Value(pol, gm.vars, F!.gmeltvars);
  if not res in F then
    return fail;
  fi;
  return res;
end);


# miscellaneous utilities
InstallMethod(PrimeField, ["IsStandardFiniteField"], 
  F -> FF(Characteristic(F), 1));


##  <#GAPDoc Label="IsFF">
##  <ManSection>
##  <Heading>Filters for standard fields</Heading>
##  <Prop Name="IsStandardPrimeField" Arg="F" />
##  <Prop Name="IsStandardFiniteField" Arg="F" />
##  <Filt Name="IsStandardFiniteFieldElement" Arg="x" Type="Category"/>
##  <Returns><K>true</K> or <K>false</K> </Returns>
##  <Description>
##  These  properties  identify  the   finite  fields  constructed  with  <Ref
##  Func="FF"  />. Prime  fields constructed  as  <C>FF(p, 1)</C>  have    the
##  property  <Ref  Prop="IsStandardPrimeField"/>.  They  are  identical  with
##  <C>GF(p)</C>,  but  calling them  via  <Ref  Func="FF"  /> we  store  some
##  additional information in these objects.
##  <P/> 
##  The fields constructed  by     <C>FF(p,k)</C> with <C>k</C> <M> &gt; 1</M> 
##  have the property <Ref Prop="IsStandardFiniteField"/>. Elements <A>x</A> in
##  such a field are in <Ref Filt="IsStandardFiniteFieldElement"/>.
##  <P/>
##  <Example>gap> F := FF(19,1);
##  GF(19)
##  gap> IsStandardFiniteField(F);
##  false
##  gap> IsStandardPrimeField(F);
##  true
##  gap> F := FF(23,48);
##  FF(23, 48)
##  gap> IsStandardFiniteField(F);
##  true
##  gap> IsStandardFiniteFieldElement(Random(F));
##  true
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallOtherMethod(IsStandardPrimeField, ["IsField"], function(K)
  local p;
  if not IsFinite(K) then
    return false;
  fi;
  p := Size(K);
  if not IsPrimeInt(p) then
    return false;
  fi;
  if K = FF(p,1) then
    return true;
  fi;
  return false;
end);

##  <#GAPDoc Label="FFElmConv">
##  <ManSection>
##  <Heading>Maps for elements of standard finite fields</Heading>
##  <Meth Name="AsVector" Label="for elements in standard finite fields"
##  Arg="a"/>
##  <Returns>a vector over prime field of <C>F</C> </Returns>
##  <Meth Name="ElementVector" Arg="F, v"/>
##  <Returns>an element in <C>F</C> </Returns>
##  <Meth Name="AsPolynomial" Label="for elements in standard finite fields"
##  Arg="a"/>
##  <Returns>a polynomial in variables of the tower of <C>F</C> </Returns>
##  <Meth Name="ElementPolynomial" Arg="F, pol"/>
##  <Returns>an element in <C>F</C> </Returns>
##  <Meth Name="SteinitzNumber" Arg="a"/>
##  <Returns>an integer</Returns>
##  <Meth Name="ElementSteinitzNumber" Arg="F, nr"/>
##  <Returns>an element in <C>F</C></Returns>
##  <Description>
##  Here,    <A>F</A>   is    always   a    standard   finite    field   (<Ref
##  Filt="IsStandardFiniteField"/>) and <A>a</A>  is an  element of  <A>F</A>.
##  <P/> 
##  <Ref Meth="AsVector" Label="for elements  in standard finite fields"/>
##  returns the coefficient  vector of <A>a</A> with respect  to the tower
##  basis of <A>F</A>. And  vice versa <Ref Meth="ElementVector"/> returns
##  the element of <A>F</A> with the given coefficient vector.
##  <P/>
##  Similarly,  <Ref Meth="AsPolynomial"  Label="for elements  in standard
##  finite   fields"/>   returns   the   (reduced)   polynomial   in   the
##  indeterminates defining  the tower of  <A>F</A>. Here, for  each prime
##  <M>r</M> dividing the degree of  the field the polynomial defining the
##  <M>k</M>-th  extension of  degree  <M>r</M> over  the  prime field  is
##  written  in the  variable  <C>x</C><M>r</M><C>_</C><M>k</M>. And  <Ref
##  Meth="ElementPolynomial"/> returns the element of <A>F</A> represented
##  by the given polynomial (which does not need to be reduced).
##  <P/>
##  Finally, <Ref  Meth="SteinitzNumber"/> returns the Steinitz  number of
##  <A>a</A>. And <Ref  Meth="ElementSteinitzNumber"/> returns the element
##  with given Steinitz number.
##  <Example><![CDATA[gap> F := FF(17, 12);
##  FF(17, 12)
##  gap> a := PrimitiveElement(F);; a := a^11-3*a^5+a;
##  ZZ(17,12,[0,1,0,0,0,14,0,0,0,0,0,1])
##  gap> v := AsVector(a);
##  < immutable compressed vector length 12 over GF(17) >
##  gap> a = ElementVector(F, v);
##  true
##  gap> ExtRepOfObj(a) = v * TowerBasis(F);
##  true
##  gap> pol := AsPolynomial(a);;
##  gap> ElementPolynomial(F, pol^10) = a^10;
##  true
##  gap> nr := SteinitzNumber(a);
##  506020624175737
##  gap> a = ElementSteinitzNumber(F, nr);
##  true
##  gap> ## primitive element of FF(17, 6)
##  gap> y := ElementSteinitzNumber(F, 17^5);
##  ZZ(17,12,[0,0,1,0,0,0,12,0,0,0,5,0])
##  gap> y = ValuePol([0,0,1,0,0,0,12,0,0,0,5,0], PrimitiveElement(F));
##  true
##  gap> x6 := Indeterminate(FF(17,1), "x6");;
##  gap> MinimalPolynomial(FF(17,1), y, x6) = DefiningPolynomial(FF(17,6));
##  true
##  ]]></Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallMethod(ElementVector, ["IsStandardFiniteField", "IsRowVector"],
function(F, vec)
  local fam, l;
  fam := FamilyObj(PrimitiveElement(F));
  l := vec * TowerBasis(F);
  if ForAll([2..Length(l)], i-> IsZero(l[i])) then
    l := l[1];
  fi;
  return ObjByExtRep(fam, l);
end);

InstallOtherMethod(SteinitzNumber, ["IsFFE"], IntFFE);
InstallOtherMethod(SteinitzNumber, 
               ["IsStandardPrimeField", "IsRingElement"],
function(F, c)
  return IntFFE(c);
end);
InstallMethod(SteinitzNumber, ["IsStandardFiniteField", "IsRingElement"],
function(F, c)
  local v;
  v := List(ExtRepOfObj(c) * PrimitivePowersInTowerBasis(F), IntFFE);
  return ValuePol(v, Characteristic(F));
end);
InstallOtherMethod(SteinitzNumber, ["IsStandardFiniteFieldElement"],
function(x)
  return SteinitzNumber(FamilyObj(x)!.wholeField, x);
end);


# reverse
InstallOtherMethod(ElementSteinitzNumber, 
                   ["IsStandardPrimeField", "IsInt"],
function(F, nr)
  return nr*One(F);
end);
InstallMethod(ElementSteinitzNumber, ["IsStandardFiniteField", "IsInt"],
function(F, nr)
  local p, fam, tc, l;
  if nr < 0 or nr >= Size(F) then
    return fail;
  fi;
  if nr = 0 then
    return Zero(F);
  fi;
  p := Characteristic(F);
  fam := FamilyObj(PrimitiveElement(F));
  tc := One(GF(p))*CoefficientsQadic(nr, p);
  ConvertToVectorRep(tc, p);
  l := tc * TowerBasis(F);
  if ForAll([2..Length(l)], i-> IsZero(l[i])) then
    l := l[1];
  fi;
  return ObjByExtRep(fam, l);
end);

InstallMethod(ViewString, ["IsStandardFiniteField"], function(F)
  return Concatenation("FF(", String(Characteristic(F)),
                       ", ", String(DegreeOverPrimeField(F)), ")");
end);
InstallMethod(ViewObj, ["IsStandardFiniteField"], function(F)
  Print(ViewString(F));
end);
InstallMethod(PrintString, ["IsStandardFiniteField"], ViewString);

InstallMethod(PrintObj, ["IsStandardFiniteField"], function(F)
  Print(PrintString(F));
end);

# helper: describe monomials in tower basis and their degrees
InstallGlobalFunction(StdMon, function(n)
  local rs, res, deg, r, k, new, ndeg, i, j;
  rs := Factors(n);
  res := [[]];
  deg := [1];
  for i in [1..Length(rs)] do
    r := rs[i];
    if i = 1 or rs[i-1] <> r then
      k := 1;
    else
      k := k+1;
    fi;
    new := ShallowCopy(res);
    ndeg := ShallowCopy(deg);
    for j in [1..r-1] do
      Append(new, List(res, a-> Concatenation(a, [r,k,j])));
      Append(ndeg, List(deg, a-> LcmInt(a,r^k)));
    od;
    res := new;
    deg := ndeg;
  od;
  return [res, deg];
end);

# just the degrees
InstallGlobalFunction(StdMonDegs, function(n)
  local f, a, res, new, i;
  if n = 1 then return [1]; fi;
  f := Collected(Factors(n));
  a := f[Length(f)];
  res := StdMonDegs(n/a[1]);
  new := List(res, l-> LcmInt(l, a[1]^a[2]));
  for i in [1..a[1]-1] do
    Append(res, new);
  od;
  return res;
end);
# map of monomials for degree n into monomials of degree m (by positions)
InstallGlobalFunction(StdMonMap, function(n, m)
  local d;
  d := StdMonDegs(m);
  return Filtered([1..Length(d)], i-> n mod d[i] = 0);
end);

# Embedding of element in FF(p,n) into FF(p,m) by Steinitz numbers
InstallGlobalFunction(EmbedSteinitz, function(p, n, m, nr)
  local l, map, c;
  if n=m or nr = 0 then return nr; fi;
  l := CoefficientsQadic(nr, p);
  map := StdMonMap(n, m){[1..Length(l)]};
  c := 0*[1..map[Length(map)]];
  c{map} := l;
  return ValuePol(c, p);
end);

InstallMethod(ExtRepOfObj,"baseFieldElm",true,
  [IsAlgebraicElement and IsAlgBFRep],0,
function(e)
  local f,l;
  f:=FamilyObj(e);
  l := 0*ShallowCopy(f!.polCoeffs{[1..f!.deg]});
  l[1] := e![1];
  MakeImmutable(l);
  return l;
end);

InstallMethod(Embedding, ["IsStandardPrimeField", "IsStandardPrimeField"],
function(K, L)
  local emb;
  if K <> L then
    Error("Embedding need identical prime fields.");
  fi;
  emb := MappingByFunction(K, L, IdFunc, false, IdFunc);
  SetFilterObj(emb, IsStandardFiniteFieldEmbedding);
  return emb;
end);
InstallMethod(Embedding, ["IsStandardPrimeField", "IsStandardFiniteField"],
function(K, L)
  local p, famL, fun, pre, emb;
  p := Characteristic(K);
  famL := FamilyObj(PrimitiveElement(L));
  fun := function(x)
    return ObjByExtRep(famL, x);
  end;
  pre := function(y)
    if IsAlgBFRep(y) then
      return y![1];
    else
      return fail;
    fi;
  end;
  emb := MappingByFunction(K, L, fun, false, pre);
  SetFilterObj(emb, IsStandardFiniteFieldEmbedding);
  return emb;
end);

##  <#GAPDoc Label="FFEmbedding">
##  <ManSection>
##  <Meth Name="Embedding" Label="for standard finite fields" Arg="H, F"/>
##  <Returns>a field homomorphism </Returns>
##  <Description>
##  Let  <A>F</A> and  <A>H</A>  be  standard finite  fields  and <A>H</A>  be
##  isomorphic to a subfield of  <A>F</A>. This function returns the canonical
##  embedding of <A>H</A> into <A>F</A>.
##  <P/> 
##  <Example>gap> F := FF(7, 360);
##  FF(7, 360)
##  gap> H := FF(7, 45);
##  FF(7, 45)
##  gap> emb := Embedding(H, F);
##  MappingByFunction( FF(7, 45), FF(7, 360), function( x ) ... end )
##  gap> y := PrimitiveElement(H);
##  ZZ(7,45,[0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
##  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
##  gap> x := y^emb;;
##  gap> ((y+One(H))^12345)^emb = (x+One(F))^12345;
##  true
##  gap> PreImageElm(emb, x^5);
##  ZZ(7,45,[0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
##  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])
##  gap> PreImageElm(emb, PrimitiveElement(F));
##  fail
##  gap> SteinitzNumber(y);
##  13841287201
##  gap> SteinitzNumber(x) mod 10^50;
##  72890819326613454654477690085519113574118965817601
##  gap> SteinitzPair(x);
##  [ 45, 13841287201 ]
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallMethod(Embedding, ["IsStandardFiniteField", "IsStandardFiniteField"],
function(K, L)
  local p, d, a, monK, monL, map, z, famL, famK, fun, imgen, fun2, pre, emb;
  p := Characteristic(K);
  if Characteristic(L) <> p then
    Error("different characteristic");
  fi;
  d := DegreeOverPrimeField(L); 
  a := DegreeOverPrimeField(K);
  if not d mod a = 0 then
    Error("not a subfield");
  fi;
  if not IsBound(K!.embcache) then
    K!.embcache := rec();
  fi;
  if IsBound(K!.embcache.(d)) then
    return K!.embcache.(d);
  fi;
  monK := StdMon(a);
  monL := StdMon(d);
  map := List(monK[1], x-> Position(monL[1], x));
  z := AsVector(Zero(L));
  famL := FamilyObj(PrimitiveElement(L));
  famK := FamilyObj(PrimitiveElement(K));
  fun := function(x)
    local v, res, l;
    v := AsVector(x);
    res := ShallowCopy(z);
    res{map} := v;
    l := res*TowerBasis(L);
    if ForAll([2..Length(l)], i-> IsZero(l[i])) then
      l := l[1];
    fi;
    return ObjByExtRep(famL, l);
  end;
##  # Maybe this version of "fun" can be made faster?  
##    imgen := fun(PrimitiveElement(K));
##    fun2 := function(x)
##      local c;
##      c := x![1];
##      if not IsList(c) then
##        return c*One(L);
##      fi;
##      return ValuePol(c, imgen);
##    end;
  pre := function(y)
    local v, l;
    v := AsVector(L, y);
    if ForAny([1..d], i-> not (i in map or IsZero(v[i]))) then
      return fail;
    fi;
    l := v{map};
    if ForAll([2..a], i-> IsZero(l[i])) then
      return ObjByExtRep(famK, l[1]);
    else
      return ObjByExtRep(famK, l*TowerBasis(K));
    fi;
  end;
  emb := MappingByFunction(K, L, fun, false, pre);
  SetFilterObj(emb, IsStandardFiniteFieldEmbedding);
  K!.embcache.(d) := emb;
  return emb;
end);
InstallMethod(IsSurjective, ["IsStandardFiniteFieldEmbedding"],
function(emb)
  return DegreeOverPrimeField(Source(emb)) = DegreeOverPrimeField(Range(emb));
end);
InstallOtherMethod(PreImageElm, 
     ["IsStandardFiniteFieldEmbedding", "IsRingElement"],
function(emb, y)
  return emb!.prefun(y);
end);

##  <#GAPDoc Label="SteinitzPair">
##  <ManSection>
##  <Oper Name="SteinitzPair" Arg="a"/>
##  <Returns>a pair of integers</Returns>
##  <Meth Name="SteinitzPair" Label="for Steinitz number" Arg="K, snr"/>
##  <Returns>a pair of integers</Returns>
##  <Meth Name="SteinitzNumber" Label="for Steinitz pair" Arg="K, pair"/>
##  <Returns>an integer</Returns>
##  <Description>
##  The    argument    <A>a</A>    must    be    an    element    in    <Ref
##  Filt="IsStandardFiniteFieldElement"/>.  Then <Ref  Oper="SteinitzPair"/>
##  returns a pair  <C>[d, nr]</C> where <C>d</C> is the  degree of <A>a</A>
##  over  the prime  field <C>FF(p,  1)</C>  and <C>nr</C>  is the  Steinitz
##  number of <A>a</A> considered as element of <C>FF(p, d)</C>.
##  <P/>
##  In the second variant a standard  finite field <A>K</A> is given and the
##  Steinitz number of an element in <A>K</A> and the result is the Steinitz
##  pair of the corresponding element.
##  <P/>
##  The inverse map  is provided by a method  for <Ref Meth="SteinitzNumber"
##  Label="for Steinitz  pair"/> which  gets a standard  finite field  and a
##  Steinitz pair.
##  <Example>gap> F := FF(7, 360);
##  FF(7, 360)
##  gap> t := ElementSteinitzNumber(F, 7^10);; # prim. elt of FF(7,12)
##  gap> sp := SteinitzPair(t);
##  [ 12, 117649 ]
##  gap> H := FF(7, 12);
##  FF(7, 12)
##  gap> b := ElementSteinitzNumber(H, 117649);
##  ZZ(7,12,[0,1,0,0,0,0,0,0,0,0,0,0])
##  gap> Value(MinimalPolynomial(FF(7,1), t), b);
##  ZZ(7,12,[0])
##  gap> nr := SteinitzNumber(t);
##  282475249
##  gap> nr = SteinitzNumber(F, sp);
##  true
##  gap> sp = SteinitzPair(F, nr);
##  true
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
# SteinitzPair(x):  [deg, nr] where deg is the degree of x over prime field 
# and nr the Steinitz number of x as element in FF(p, deg).
BindGlobal("StPairIntVec", function(K, v)
  local p, d, map, ind, sd, vv;
  p := Characteristic(K);
  d := DegreeOverPrimeField(K);
  map := StdMonDegs(d);
  ind := Filtered([1..d], i-> IsBound(v[i]) and not IsZero(v[i]));
  if Length(ind) = 0 then
    return [1,0];
  fi;
  sd := Lcm(Set(map{ind}));
  ind := Filtered([1..d], i-> sd mod map[i] = 0);
  while ind[Length(ind)] > Length(v) do
    Remove(ind);
  od;
  vv := v{ind};
  return [sd, ValuePol(vv, p)];
end);
GAPInfo.tmpmeth := function(F, x)
  return StPairIntVec(F, List(AsVector(x), IntFFE));
end;
InstallOtherMethod(SteinitzPair, 
    ["IsStandardFiniteField", "IsStandardFiniteFieldElement"], GAPInfo.tmpmeth);
Unbind(GAPInfo.tmpmeth);
InstallOtherMethod(SteinitzPair, ["IsStandardPrimeField", "IsFFE"], 
function(F, x)
  return [1, SteinitzNumber(F,x)];
end);

# with Steinitz number instead of element
GAPInfo.tmpmeth := function(F, nr)
  local p;
  p := Characteristic(F);
  return StPairIntVec(F, CoefficientsQadic(nr, p));
end;
InstallOtherMethod(SteinitzPair, 
    ["IsStandardFiniteField", "IsInt"], GAPInfo.tmpmeth);
Unbind(GAPInfo.tmpmeth);
InstallOtherMethod(SteinitzPair, ["IsStandardPrimeField", "IsInt"], 
function(F, nr)
  return [1, nr];
end);

InstallMethod(SteinitzPair, ["IsStandardFiniteFieldElement"], 
function(x)
  return SteinitzPair(FamilyObj(x)!.wholeField, x);
end);
# works only for prime fields
InstallOtherMethod(SteinitzPair, ["IsStandardPrimeField","IsFFE"], 
{F, x}-> [1, IntFFE(x)]); 
InstallOtherMethod(SteinitzPair, ["IsFFE"], x-> [1, IntFFE(x)]); 

# convert Steinitz pair in Steinitz number
GAPInfo.tmp := 
function(K, st)
  local p, d, v, map, ind, vv;
  p := Characteristic(K);
  d := DegreeOverPrimeField(K);
  if st[1] = d then
    return st[2];
  fi;
  v := CoefficientsQadic(st[2], p);
  map := StdMonDegs(d);
  # indices of subsequence of tower basis of subfield
  ind := Filtered([1..d], i-> st[1] mod map[i] = 0);
  if Length(ind) > Length(v) then
    ind := ind{[1..Length(v)]};
  fi;
  vv := 0*[1..d];
  vv{ind} := v;
  return ValuePol(vv, p);
end;
InstallOtherMethod(SteinitzNumber, ["IsStandardFiniteField", "IsList"],
GAPInfo.tmp);
InstallOtherMethod(SteinitzNumber, ["IsStandardPrimeField", "IsList"],
function(K, pair)
  if pair[1] <> 1 then 
    Error("Steinitz pair not in prime field.");
  fi;
  return pair[2];
end);

##  <#GAPDoc Label="SteinitzPairConwayGenerator">
##  <ManSection>
##  <Func Name="SteinitzPairConwayGenerator" Arg="F"/>
##  <Returns>a pair of integers </Returns>
##  <Description>
##  For a standard finite field <A>F</A>  of order <M>q</M> for which a Conway
##  polynomial  (see <Ref  BookName="Reference" Func="ConwayPolynomial"/>)  is
##  known this function returns the <Ref Oper="SteinitzPair"/> for the element
##  of <A>F</A> corresponding to <C>Z(q)</C>  (which is by definition the zero
##  of  the Conway  polynomial in  <A>F</A> with  the smallest Steinitz number
##  which is compatible with the choice in all proper subfields).
##  <P/> 
##  This is used to construct the
##  <Ref Func="StandardIsomorphismGF"/> for <A>F</A>.
##  <Example>gap> F := FF(23,18);
##  FF(23, 18)
##  gap> st := SteinitzPairConwayGenerator(F);
##  [ 18, 1362020736983803830549380 ]
##  gap> st9 := SteinitzPairConwayGenerator(FF(23,9));
##  [ 9, 206098743447 ]
##  gap> st6 := SteinitzPairConwayGenerator(FF(23,6));
##  [ 6, 45400540 ]
##  gap> z  := ElementSteinitzNumber(F, st[2]);;
##  gap> z9 := ElementSteinitzNumber(F, SteinitzNumber(F, st9));;
##  gap> z6 := ElementSteinitzNumber(F, SteinitzNumber(F, st6));;
##  gap> e9 := (Size(F)-1)/(23^9-1);
##  1801152661464
##  gap> e6 := (Size(F)-1)/(23^6-1);
##  21914624580056211
##  gap> z9 = z^e9;
##  true
##  gap> z6 = z^e6;
##  true
##  gap> l := Filtered(ZeroesConway(F), x-> x^e9 = z9 and x^e6 = z6);;
##  gap> List(l, SteinitzNumber);
##  [ 1362020736983803830549380 ]
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallOtherMethod(SteinitzPairConwayGenerator, ["IsStandardPrimeField"],
function(K)
  local pl;
  pl := ConwayPol(Characteristic(K), 1) * One(K);
  return [1, SteinitzNumber(K, -pl[1])];
end);

InstallMethod(SteinitzPairConwayGenerator, ["IsStandardFiniteField"],
function(K)
  local d, p, rs, compat, exp, L, c, emb, pl, z, st, i, r;
  d := DegreeOverPrimeField(K);
  p := Characteristic(K);
  if not IsCheapConwayPolynomial(p, d) then
    return fail;
  fi;
  rs := Set(Factors(d));
  compat := [];
  exp := [];
  for r in rs do
    L := FF(p, d/r);
    c := SteinitzPairConwayGenerator(L);
    emb := Embedding(L, K);
    Add(compat, Image(emb, ElementSteinitzNumber(L,c[2])));
    Add(exp, (Size(K)-1)/(Size(L)-1));
  od;

  pl := ConwayPol(p, d) * One(K);
  #z := FindConjugateZeroes(K, pl, p);
  z := ZeroesConway(K);
  st := List(z, c-> SteinitzNumber(K, c));
  SortParallel(st, z);
  # now find Steinitz smallest zero compatible with compat
  i := First([1..d], i-> ForAll([1..Length(rs)], j-> z[i]^exp[j] = compat[j]));
  return SteinitzPair(K, z[i]);
end);

##  <#GAPDoc Label="StandardIsomorphismGF">
##  <ManSection>
##  <Func Name="StandardIsomorphismGF" Arg="F" />
##  <Returns>a field isomorphism </Returns>
##  <Description>
##  The   argument  <A>F</A>   must  be   a  standard   finite  field,   say
##  <C>FF(p,n)</C>  such  that  &GAP;   can  generate  <C>GF(p,n)</C>.  This
##  function returns the field  isomorphism from <C>GF(p,n)</C> to <A>F</A>,
##  which sends <C>Z(p,n)</C> to the  element with Steinitz pair computed by
##  <Ref Func="SteinitzPairConwayGenerator" />.
##  <P/> 
##  <Example>gap> F := FF(13,21);
##  FF(13, 21)
##  gap> iso := StandardIsomorphismGF(F);
##  MappingByFunction( GF(13^21), FF(13, 21), function( x ) ... end )
##  gap> K := GF(13,21);
##  GF(13^21)
##  gap> x := Random(K);;
##  gap> l := [1,2,3,4,5];;
##  gap> ValuePol(l, x)^iso = ValuePol(l, x^iso);
##  true
##  gap> y :=  ElementSteinitzNumber(F, SteinitzPairConwayGenerator(F)[2]);;
##  gap> PreImageElm(iso, y);
##  z
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallGlobalFunction(StandardIsomorphismGF, 
function(L)
  local p, d, K, z, fun, a, mat, imat, zz, pre, emb, i;
  p := Characteristic(L);
  d := DegreeOverPrimeField(L);
  if IsBound(L!.isocache) then
    return L!.isocache;
  fi;
  if d = 1 then
    emb := IdentityMapping(L);
    SetFilterObj(emb, IsFieldHomomorphism and IsBijective);
    return emb;
  fi;
  K := GF(p,d);
  z := ElementSteinitzNumber(L, SteinitzPairConwayGenerator(L)[2]);
  fun := function(x)
    local c, bas;
    if IsPositionalObjectRep(x) then
      c := x![1];
    else
      bas := CanonicalBasis(K);
      c := Coefficients(bas, x);
    fi;
    return ValuePol(c, z);
  end;
  a := One(L);
  mat := [];
  for i in [1..d] do
    Add(mat, ExtRepOfObj(a));
    a := a*z;
  od;
  imat := mat^-1;
  zz := Z(Size(K));
  pre := function(y)
    return ValuePol(ExtRepOfObj(y)*imat, zz);
  end;
  emb := MappingByFunction(K, L, fun, false, pre);
  SetFilterObj(emb, IsFieldHomomorphism and IsBijective);
  L!.isocache := emb;
  return emb;
end);

##  <#GAPDoc Label="FFArith">
##  <ManSection>
##  <Oper Name="ZZ" Arg="p, n, coeffs" />
##  <Oper Name="ZZ" Arg="p, n, ffe" Label="for IsFFE" />
##  <Returns>an element in <C>FF(<A>p</A>, <A>n</A>)</C></Returns>
##  <Description>
##  For a  prime <A>p</A>, positive  integer <A>n</A> and an  integer list
##  <A>coeffs</A>  this function  returns the  element in  <C>FF(<A>p</A>,
##  <A>n</A>)</C>  represented by  the  polynomial  with coefficient  list
##  <A>coeffs</A> modulo <A>p</A>. Elements  in standard finite fields are
##  also printed in this way.
##  <P/>
##  For convenience the third argument <A>ffe</A> can be in `GF(p,n)` (see
##  <Ref  BookName="Reference"  Func="GF"  Label="for  characteristic  and
##  degree" /> and <Ref BookName="Reference" Filt="IsFFE"/>). This returns
##  the image of <A>ffe</A>  under the <Ref Func="StandardIsomorphismGF"/>
##  of <C>FF(<A>p</A>,<A>n</A>)</C>.
##  <Example>gap> x := ZZ(19,5,[1,2,3,4,5]);
##  ZZ(19,5,[1,2,3,4,5])
##  gap> a := PrimitiveElement(FF(19,5));
##  ZZ(19,5,[0,1,0,0,0])
##  gap> x = [1,2,3,4,5]*[a^0,a^1,a^2,a^3,a^4];
##  true
##  gap> One(FF(19,5)); # elements in prime field abbreviated
##  ZZ(19,5,[1])
##  gap> One(FF(19,5)) = ZZ(19,5,[1]);
##  true
##  gap> ZZ(19,5,Z(19^5)); # zero of ConwayPolynomial(19,5)
##  ZZ(19,5,[12,5,3,4,5])
##  </Example>
##  </Description>
##  </ManSection>
##  
##  <ManSection>
##  <Func Name="MoveToSmallestStandardField" Arg="x" />
##  <Meth Name="\+" Arg="x, y" Label="for standard finite field elements" />
##  <Meth Name="\*" Arg="x, y" Label="for standard finite field elements" />
##  <Meth Name="\-" Arg="x, y" Label="for standard finite field elements" />
##  <Meth Name="\/" Arg="x, y" Label="for standard finite field elements" />
##  <Returns>a field element </Returns>
##  <Description>
##  Here <A>x</A> and <A>y</A> must  be elements in standard finite fields
##  (of the same characteristic).
##  <P/> 
##  Then  <Ref Func="MoveToSmallestStandardField"  /> returns  the element
##  <A>x</A> as element of the smallest possible degree extension over the
##  prime field.
##  <P/>
##  The arithmetic operations are even possible when <A>x</A> and <A>y</A>
##  are not  represented as elements in  the same field. In  this case the
##  elements are first mapped to the smallest field containing both.
##  <Example>gap> F := FF(1009,4);
##  FF(1009, 4)
##  gap> G := FF(1009,6);
##  FF(1009, 6)
##  gap> x := (PrimitiveElement(F)+One(F))^13;
##  ZZ(1009,4,[556,124,281,122])
##  gap> y := (PrimitiveElement(G)+One(G))^5;
##  ZZ(1009,6,[1,5,10,10,5,1])
##  gap> x+y;
##  ZZ(1009,12,[557,0,936,713,332,0,462,0,843,191,797,0])
##  gap> x-y;
##  ZZ(1009,12,[555,0,73,713,677,0,97,0,166,191,212,0])
##  gap> x*y;
##  ZZ(1009,12,[253,289,700,311,109,851,345,408,813,657,147,887])
##  gap> x/y;
##  ZZ(1009,12,[690,599,714,648,184,217,563,130,251,675,73,782])
##  gap> z  := -y + (x+y);
##  ZZ(1009,12,[556,0,0,713,0,0,784,0,0,191,0,0])
##  gap> SteinitzPair(z);
##  [ 4, 125450261067 ]
##  gap> x=z;
##  true
##  gap> MoveToSmallestStandardField(z);
##  ZZ(1009,4,[556,124,281,122])
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
# creating finite field elements
InstallMethod(ZZ, ["IsPosInt", "IsPosInt", "IsList"], 
function(p,d,c)
  local F, i, fam;
  F := FF(p,d);
  if ForAll([2..Length(c)], i-> c[i] = 0) then
    return c[1] * One(F);
  fi;
  if Length(c) < d then
    c := ShallowCopy(c);
    for i in [Length(c)+1..d] do
      c[i] := 0;
    od;
  fi;
  fam := FamilyObj(PrimitiveElement(F));
  c := [One(FF(p,1)) * c];
  ConvertToVectorRep(c[1], p);
  Objectify(fam!.extType, c);
  return c;
end);
# utility method to embed FFE into standard field
InstallOtherMethod(ZZ, ["IsPosInt", "IsPosInt", "IsFFE"],
function(p,d,x)
  local F, emb;
  F := FF(p,d);
  emb := StandardIsomorphismGF(F);
  return x^emb;
end);

# nicer print/view for standard finite field elements
InstallMethod(PrintString, ["IsStandardFiniteFieldElement"],
function(x)
  local F, c, res;
  F := FamilyObj(x)!.wholeField;
  c := x![1];
  if IsFFE(c) then c := [c]; fi;
  res := Concatenation( "ZZ(", String(Characteristic(F)),",",
            String(DegreeOverPrimeField(F)),",",String(List(c, IntFFE)),")");
  RemoveCharacters(res, " ");
  res := SubstitutionSublist(res, ",", "\<,\>");
  return res;
end);
InstallMethod(PrintObj, ["IsStandardFiniteFieldElement"], 
              function(x) Print(PrintString(x)); end);
InstallMethod(ViewString, ["IsStandardFiniteFieldElement"], PrintString);
InstallMethod(String, ["IsStandardFiniteFieldElement"], 
              x-> StripLineBreakCharacters(PrintString(x)));

# arithmetic for elements from different standard fields
InstallGlobalFunction(MoveToCommonStandardField, function(a, b)
  local Fa, Fb, p, da, db, d, F;
  Fa := FamilyObj(a)!.wholeField;
  Fb := FamilyObj(b)!.wholeField;
  p := Characteristic(Fa);
  if p <> Characteristic(Fb) then Error("different characteristic"); fi;
  da := DegreeOverPrimeField(Fa);
  db := DegreeOverPrimeField(Fb);
  d := LcmInt(da, db);
  F := FF(p, d);
  if da < d then
    a := a^Embedding(Fa, F);
  fi;
  if db < d then
    b := b^Embedding(Fb, F);
  fi;
  return [a, b];
end);
InstallMethod(\+, IsNotIdenticalObj, ["IsStandardFiniteFieldElement",
"IsStandardFiniteFieldElement"], function(a, b)
  return CallFuncList(\+, MoveToCommonStandardField(a, b));
end);
InstallMethod(\*, IsNotIdenticalObj, ["IsStandardFiniteFieldElement",
"IsStandardFiniteFieldElement"], function(a, b)
  return CallFuncList(\*, MoveToCommonStandardField(a, b));
end);
InstallMethod(\=, IsNotIdenticalObj, ["IsStandardFiniteFieldElement",
"IsStandardFiniteFieldElement"], function(a, b)
  return Characteristic(a) = Characteristic(b) and
         SteinitzPair(a) = SteinitzPair(b);
end);
# utility to represent element in smallest field
InstallGlobalFunction(MoveToSmallestStandardField, function(a)
  local Fa, d, st, F;
  if IsFFE(a) then
    # also handle elements in GF(p,n)
    d := DegreeOverPrimeField(DefaultField(a));
    if d = 1 then return a; fi;
    a := a^StandardIsomorphismGF(FF(Characteristic(a), d));
  fi;
  Fa := FamilyObj(a)!.wholeField;
  d := DegreeOverPrimeField(Fa);
  st := SteinitzPair(a);
  if st[1] = d then
    return a;
  else
    F := FF(Characteristic(Fa), st[1]);
    return PreImageElm(Embedding(F, Fa), a);
  fi;
end);

