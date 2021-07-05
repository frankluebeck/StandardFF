#############################################################################
##
#A  StandardFF.gi                                                Frank Lübeck
##  
##  The files StandardFF.g{i,d}  contain code to compute  standard models of
##  finite fields.  Fields are  returned as field  towers with  prime degree
##  indices, and as single extensions over the prime field.
##  

# we cache known extensions in prime field
SFFHelper.PrepareGFpForStandardFF := function(p, r, variant...)
  local nam, F, l, lr;
  if "Trace" in variant then
    nam := "rPowerExtTrace";
  else
    nam := "rPowerExt";
  fi;

  F := GF(p);
  SetIsStandardPrimeField(F, true);
  if not IsBound(F!.(nam)) then
    l := [];
    F!.(nam) := l;
  else
    l := F!.(nam);
  fi;
  if not IsBound(l[r]) then
    lr := [];
    l[r] := lr;
  else
    lr := l[r];
  fi;
  return lr;
end;

# k-th iterated polynomials of degree p
SFFHelper.ArtinSchreierStandardFF := function(p, k)
  local F, lr, v;
  F := GF(p);
  SetIsStandardPrimeField(F, true);
  lr := SFFHelper.PrepareGFpForStandardFF(p, p);
  if IsBound(lr[k]) then 
    return lr[k];
  fi;
  if k=1 then
    # X^P - X -1 is irreducible
    v := Indeterminate(F, Concatenation("x",String(p),"_1"));
    Add(lr, [v^p-v-1,v]);
  else
    v := Indeterminate(F, Concatenation("x",String(p),"_",String(k)));
    SFFHelper.ArtinSchreierStandardFF(p, k-1);
    Add(lr, [v^p - v - Product(List(lr, a-> a[2]))^(p-1), v]);
  fi;
  return lr[k];
end;

# for odd p extensions of 2 power degree
SFFHelper.2PowerStandardFF := function(p, k)
  local F, lr, p2, v, l, m, a, b, q, q2, pol1, v1, c1, n1, c2, n2, u;
  F := GF(p);
  SetIsStandardPrimeField(F, true);
  lr := SFFHelper.PrepareGFpForStandardFF(p, 2);
  if IsBound(lr[k]) then 
    return lr[k];
  fi;
  p2 := (p-1)/2;
  v := Indeterminate(F, Concatenation("x2_", String(k)));
  if k = 1 then
    v := Indeterminate(F, "x2_1");
    if p2 mod 2 = 1 then
      Add(lr, [v^2+1, v]);
    else
      l := 2;
      p2 := p2/2;
      while p2 mod 2 = 0 do
        l := l+1;
        p2 := p2/2;
      od;
      # find random 2^l-th root of 1, half of elements will do
      m := One(F);
      while IsZero(m) or IsOne(m) do
        a := Random(F);
        m := a^((p-1)/2);
      od;
      a := a^((p-1)/2^l);
      # now change a for u=l-2 .. 0 such that a^(2^u) is the smaller root
      # of a^(2^(u+1))
      for u in [l-2,l-3..0] do
        b := a^(2^u);
        if IntFFE(-b) < IntFFE(b) then
          a := a^(1+2^(l-1-u));
        fi;
      od;
      Add(lr, [v^2-a, v]);
    fi;
  elif k=2 then
    SFFHelper.2PowerStandardFF(p, 1);
    if p2 mod 2 = 1 then
      # we proceed as above but now for q = p^2 instead of p
      q := p^2;
      q2 := (q-1)/4;
      l := 2;
      while q2 mod 2 = 0 do
        l := l+1;
        q2 := q2/2;
      od;
      pol1 := lr[1][1];
      v1 := lr[1][2];
      m := v1^0;
      while IsOne(m) or IsZero(m) do
        a := Random(F)+Random(F)*v1;
        m := PowerMod(a, (q-1)/2, pol1);
      od;
      a := PowerMod(a, (q-1)/2^l, pol1);
      # here we use Steinitz number x + y*p for x+y*v1 for comparison
      for u in [l-2,l-3..0] do
        b := PowerMod(a, 2^u, pol1);
        c1 := List(CoefficientsOfUnivariatePolynomial(b), IntFFE);
        n1 := c1[1]+p*c1[2];
        b := (b * m) mod pol1;
        c2 := List(CoefficientsOfUnivariatePolynomial(b), IntFFE);
        n2 := c2[1]+p*c2[2];
        if n2 < n1 then
          a := PowerMod(a, 1+2^(l-1-u), pol1);
        fi;
      od;
      Add(lr, [v^2-a, v]);
    else
      Add(lr, [v^2-lr[1][2], v]);
    fi;
  else
    SFFHelper.2PowerStandardFF(p, k-1);
    Add(lr, [v^2-lr[k-1][2], v]);
  fi;
  return lr[k];
end;

# input: r a prime and q a prime power
# output: list of factors of r-th cyclotomic polynomial over GF(q)
SFFHelper.GFCyclotomicPolynomials := function(r, q)
  local F, m, pol, mulmod, powmod, randpol, pols, res, rp, new, g, pl;
  F := GF(q);
  SetIsStandardPrimeField(F, true);
  # in characteristic 2 we fall back to the algorithm in GAP
  if q mod 2 =0 then
    return Factors(CyclotomicPolynomial(F, r));
  fi;
  # degree of irreducible factors
  m := OrderMod(q, r);
  # we compute with coefficient lists 
  pol := List([1..r], i->Z(q)^0);
  ConvertToVectorRep(pol, q);
  if m = r-1 then
    return [UnivariatePolynomial(F, pol)];
  fi;
  # multiply and power modulo X^r - 1
  mulmod := function(v, w)
    local prd, i;
    prd := ProductCoeffs(v,w);
    for i in [r+1..Length(prd)] do
      prd[i-r] := prd[i-r] + prd[i];
    od;
    return prd{[1..Minimum(Length(prd), r)]};
  end;
  powmod := function(v, k)
    local vp, res;
    vp := ShallowCopy(v);
    res := [Z(q)^0];
    while k > 0 do
      if k mod 2 = 1 then
        res := mulmod(res, vp);
        k := (k-1)/2;
      else
        k := k/2;
      fi;
      if k > 0 then
        vp := mulmod(vp, vp);
      fi;
    od;
    return res;
  end;

  # takes random element mod x^r-1, computes the trace and
  # raises to the (q-1)/2 power; result +/- 1 is likely to have
  # non-trivial Gcd with cyclotomic polynomial (see Cantor-Zassenhaus)
  randpol := function()
    local v, n, k, i, j;
    # coeffs of random element mod (x^r - 1)
    v := List([1..r], i-> Random(F));
    # trace
    n := ShallowCopy(v);
    for i in [0..r-1] do
      k := i;
      for j in [1..m-1] do
        k := (k*q) mod r;
        n[k+1] := n[k+1]+v[i+1];
      od;
    od;
    return powmod(n, (q-1)/2);
  end;

  pols := [pol];
  res := [];
  # compute random polys as above and take Gcd with known factors
  # (factors are irreducible if of degree m)
  while Length(pols) > 0 do
    rp := randpol();
    rp[1] := rp[1]-Z(q)^0;

    new := [];
    for pl in pols do
      g := GcdCoeffs(rp, pl);
      if Length(g) > 1 and Length(g) < Length(pl) then
        if Length(g) = m+1 then
          Add(res, g);
        else
          Add(new, g);
        fi;
        g := QuotRemCoeffs(pl, g)[1];
        if Length(g) = m+1 then
          Add(res, g);
        else
          Add(new, g);
        fi;
      else
        Add(new, pl);
      fi;
    od;
    pols := new;
  od;
  return List(res, v-> UnivariatePolynomial(F, v));
end;

SFFHelper.rPowerStandardFFTrace := function(p, r, k)
  local F, lr, t, a, zeta, zi, min, minmul, rr, m, z, cps, pp, vcps, 
        i, L, zc, q, j, sn, snzi, v, cc, e, qr, vec, st, mp, c, K;
  F := GF(p);
  SetIsStandardPrimeField(F, true);
  # let z be the smallest primitive r-th root of
  # degree m over F_p
  # we store in lr[k+1] a list with 5 entries:
  #   * poly for k-th extension K of degree r
  #   * its variable 
  #   * finite field tower for K - K[z]
  #   * a in K[z]
  #   * t such that a is of order r^t and has no r-th root
  lr := SFFHelper.PrepareGFpForStandardFF(p, r, "Trace");
  if IsBound(lr[k+1]) then 
    return lr[k+1];
  fi;
  if k = 0 then
    if (p-1) mod r = 0 then
      # simplest case, prime field has r-th roots of 1
      t := 1;
      a := (p-1)/r;
      while a mod r = 0 do
        t := t+1; 
        a := a/r;
      od;
      # random primitive r^t-th root, (1-1/r) of all elements are good
      a := Random(F);
      zeta := a^((p-1)/r);
      while IsOne(zeta) or IsZero(zeta) do
        a := Random(F);
        zeta := a^((p-1)/r);
      od;
      a := a^((p-1)/r^t);
      # now find the Steinitz lexicographically smallest a
      for i in [t-1,t-2..0] do
        zi := a^(r^i);
        min := zi;
        minmul := 0;
        if i=t-1 then
          rr := r-2;
        else
          rr := r-1;
        fi;
        for j in [1..rr] do
          zi := zi*zeta;
          if IntFFE(zi) < IntFFE(min) then
            min := zi;
            minmul := j;
          fi;
        od;
        a := a^(1+minmul*r^(t-1-i));
        if i=t-1 then
          zeta := min;
        fi;
      od;
      Add(lr, [,a,F,a,t]);
    else
      # need to extend to minimal field with r-th roots of 1
      m := OrderMod(p, r);
      z := Indeterminate(F, "z");
      cps := List(SFFHelper.GFCyclotomicPolynomials(r, p), pl-> Value(pl, z));
      pp := [1];
      # find and use Steinitz smallest poly in cps
      for i in [1..m] do Add(pp, pp[i]*p); od;
      vcps := List(cps, pl-> List(CoefficientsOfUnivariatePolynomial(pl),
              IntFFE)*pp);
      i := Position(vcps, Minimum(vcps));
      cps := cps[i];
      L := FiniteFieldTower(p, cps);
      zc := GeneratorsOfField(L)[1];
      q := p^m;
      # r-valuation in q-1
      t := 1;
      a := (q-1)/r;
      while a mod r = 0 do
        t := t+1; 
        a := a/r;
      od;
      if t = 1 then
        a := zc;
      else
        # random primitive r^t-th root, (1-1/r) of all elements are good
        a := Random(L);
        zeta := a^((q-1)/r);
        while IsOne(zeta) or IsZero(zeta) do
          a := Random(L);
          zeta := a^((q-1)/r);
        od;
        a := a^((q-1)/r^t);
        # adjust such that zeta = zc
        j := 1;
        zi := zeta;
        while zi <> zc do
          j := j+1;
          zi := zi*zeta;
        od;
        if j <> 1 then
          a := a^j;
          zeta := zc;
        fi;
        # now find the Steinitz lexicographically smallest a
        for i in [t-2,t-3..0] do
          zi := a^(r^i);
          min := zi;
          minmul := 0;
          sn := SteinitzNumber(zi);
          for j in [1..r-1] do
            zi := zi*zeta;
            snzi := SteinitzNumber(zi);
            if snzi < sn then
              min := zi;
              minmul := j;
              sn := snzi;
            fi;
          od;
          a := a^(1+minmul*r^(t-1-i));
        od;
      fi;
      # we remember cps in first position for reuse
      Add(lr, [cps, a!.rep, L, a, t]);
    fi;
  else
    # recursion
    SFFHelper.rPowerStandardFFTrace(p, r, k-1);
    L := lr[k][3];
    a := lr[k][4];
    t := lr[k][5];
    v := Indeterminate(F, Concatenation("x",String(r),"_",String(k)));
    m := OrderMod(p, r);
    if m = 1 then
      # easy case with r-th roots of 1 in prime field
      lr[k+1] := [v^r-lr[k][2], v, 0, v, t+1];
    else
      # we now compute in L[X]/(X^r-a), represented as vectors/L
      # c is r-th root of a
      # the new generator cc is the trace of c in degree m subfield
      cc := (0*a)*[1..r];
      cc[2] := a^0;
      for i in [1..m-1] do
        # for q = p^(r^k) compute c^(q^i), we use c^(r^(t+1)) = 1 and c^r = a
        e := PowerMod(p,i*r^k, r^(t+1));
        qr := QuotientRemainder(e, r);
        cc[qr[2]+1] := cc[qr[2]+1] + a^qr[1];
      od;
      # compute minimal polynomial of cc over L, will have coefficients
      # in degree m subfield K of L
      # [Remark: here we compute the minimal polynomial of the vector of 1
      # Starting with a random vector needs one more multiplication,
      # but we can take random r columns (corresp. to r elements of the
      # K-basis to compute the minimal polynomial.
      vec := 0*cc;
      vec[1] := a^0;
      st := FindLinearCombination(vec, false);
      vec := cc;
      st := FindLinearCombination(vec, st[1]);
      for i in [2..r] do
        # the powers cc^j, j = 2..r
        vec := ProductCoeffs(vec, cc);
        for j in [r+1..Length(vec)] do
          vec[j-r] := vec[j-r] + vec[j]*a;
        od;
        vec := vec{[1..r]};
        st := FindLinearCombination(vec, st[1]);
      od;
      mp := v^0;
      for i in [r,r-1..1] do
        mp := mp*v-PolynomialFiniteFieldTowerElement(st[2][i]);
      od;
      # now express c a L-linear combination of powers of cc = v
      c := 0*vec;
      c[2] := a^0;
      st := FindLinearCombination(c, st[1]);
      # extension of degree r^(k-1), the new L and the new a
      K := PreviousFiniteFieldTower(L);
      # shift coefficients in st[2] to larger extension
      # was written in K[z][c], now rewritten in K[cc][z]
      a := NullMat(m, r, K);
      for i in [1..m] do 
        for j in [1..r] do
          if IsBound(st[2][j]!.rep[i]) then
            a[i][j] := st[2][j]!.rep[i];
          fi;
        od;
      od;
      # last is stored cyclotomic polynomial
      L := ExtendFiniteFieldTower(K, mp, lr[1][1]);
      K := PreviousFiniteFieldTower(L);
      a := List(a, row-> FiniteFieldTowerElement(K, row));
      a := FiniteFieldTowerElement(L, a);

      lr[k+1] := [mp, v, L, a, t+1];
    fi;
  fi;
  return lr[k+1];
end;

# we cache the standard irreducible polynomials of prime degree
# via the Steinitz number of the polynomial minus the leading term.
SFFHelper.StandardIrreducibleCoeffListCACHE := function(L, r, a)
  local p, k, res, v;
  if UserPreference("StandardFFUseCache") <> true then
    return StandardIrreducibleCoeffList(L, r, a);
  fi;
  if not IsBound(SFFHelper.IRRPRIMDEGCACHE) then
    ReadPackage("StandardFF", "data/IRRPRIMDEGCACHE");
  fi;
  p := Characteristic(L);
  k := LogInt(DegreeOverPrimeField(L), r);
  if p > 10000 then
    return StandardIrreducibleCoeffList(L, r, a);
  fi;
  if IsBound(SFFHelper.IRRPRIMDEGCACHE[p]) then
    v := First(SFFHelper.IRRPRIMDEGCACHE[p], l-> l[1]=r and l[2]=k);
    if v <> fail then
      res := CoefficientsQadic(v[3], Size(L)); 
      while Length(res) < r do
        Add(res, 0);
      od;
      Add(res, 1);
      return List(res, i-> ElementSteinitzNumber(L,i));
    fi;
  else
    SFFHelper.IRRPRIMDEGCACHE[p] := [];
  fi;
  res := StandardIrreducibleCoeffList(L, r, a);
  v := res{[1..r]};
  v := List(v, SteinitzNumber);
  v := ValuePol(v, Size(L));
  Add(SFFHelper.IRRPRIMDEGCACHE[p], [r, k, v]);
  return res;
end;

SFFHelper.rPowerStandardFFRandom := function(p, r, k)
  local F, lr, t, a, zeta, zi, min, minmul, rr, L, 
        v, m, vec, mp, K, i, j, Lst, ast;
  F := GF(p);
  SetIsStandardPrimeField(F, true);
  lr := SFFHelper.PrepareGFpForStandardFF(p, r);
  if IsBound(lr[k+1]) then 
    return lr[k+1];
  fi;
  if k = 0 then
    if (p-1) mod r = 0 then
      # simplest case, prime field has r-th roots of 1
      t := 1;
      a := (p-1)/r;
      while a mod r = 0 do
        t := t+1; 
        a := a/r;
      od;
      # random primitive r^t-th root, (1-1/r) of all elements are good
      a := Random(F);
      zeta := a^((p-1)/r);
      while IsOne(zeta) or IsZero(zeta) do
        a := Random(F);
        zeta := a^((p-1)/r);
      od;
      a := a^((p-1)/r^t);
      # now find the Steinitz lexicographically smallest a
      for i in [t-1,t-2..0] do
        zi := a^(r^i);
        min := zi;
        minmul := 0;
        if i=t-1 then
          rr := r-2;
        else
          rr := r-1;
        fi;
        for j in [1..rr] do
          zi := zi*zeta;
          if IntFFE(zi) < IntFFE(min) then
            min := zi;
            minmul := j;
          fi;
        od;
        a := a^(1+minmul*r^(t-1-i));
        if i=t-1 then
          zeta := min;
        fi;
      od;
      Add(lr, [,a,F,a,t]);
    else
      # we use smallest/random polynomials, here we just fix 1 as
      # norm of next step
      Add(lr, [, , GF(p), -One(GF(p))]);
    fi;
  else
    # recursion
    SFFHelper.rPowerStandardFFRandom(p, r, k-1);
    v := Indeterminate(F, Concatenation("x",String(r),"_",String(k)));
    if  (p-1) mod r = 0 then
      # easy case with r-th roots of 1 in prime field
      t := lr[k][5];
      lr[k+1] := [v^r-lr[k][2], v, 0, v, t+1];
    else
      L := lr[k][3];
      a := lr[k][4];
      #if IsFiniteFieldTower(L) and IsFiniteFieldTower(L!.prev) then
      if IsFiniteFieldTower(L) then
        # go to simple extension for efficiency
        Lst := StandardFiniteField(p, r^(k-1));
        ast := ElementSteinitzNumber(Lst, SteinitzNumber(a));
        #vec := StandardIrreducibleCoeffList(Lst, r, ast);
        vec := SFFHelper.StandardIrreducibleCoeffListCACHE(Lst, r, ast);
        vec := List(vec, c-> ElementSteinitzNumber(L, SteinitzNumber(Lst,c)));
      else
        #vec := StandardIrreducibleCoeffList(L, r, a);
        vec := SFFHelper.StandardIrreducibleCoeffListCACHE(L, r, a);
      fi;
      mp := v^0;
      for i in [r,r-1..1] do
        mp := mp*v+PolynomialFiniteFieldTowerElement(vec[i]);
      od;
      K := ExtendFiniteFieldTower(L, mp);
      lr[k+1] := [mp, v, K, -ElementSteinitzNumber(K, Size(L))];
    fi;
  fi;
  return lr[k+1];
end;
#MakeImmutable(SFFHelper);

##  <#GAPDoc Label="FF">
##  <ManSection>
##  <Heading>Constructing standard finite fields</Heading>
##  <Func Name="StandardFiniteField" Arg="p, n"/>
##  <Func Name="FF" Arg="p, n"/>
##  <Returns>a finite field</Returns>
##  <Func Name="StandardFiniteFieldTower" Arg="p, n"/>
##  <Returns>a tower of finite fields</Returns>
##  <Attr Name="Tower" Arg="F"/>
##  <Returns>a tower of finite fields</Returns>
##  <Func Name="StandardPrimeDegreePolynomial" Arg="p, r, k" />
##  <Returns>a polynomial of degree <A>r</A></Returns>
##  <Description>
##  The   arguments   are   a   prime  <A>p</A>   and   a   positive   integer
##  <A>n</A>.   The   function  <Ref   Func="FF"/>   (or   its  synomym   <Ref
##  Func="StandardFiniteField"/>)  is  one  of  the  main  functions  of  this
##  package.   It   returns   a   standardized   field   <C>F</C>   of   order
##  <M><A>p</A>^{<A>n</A>}</M>.  It   is  implemented  as  a   simple  extension
##  over  the   prime  field  <C>GF(p)</C>  using   <Ref  BookName="Reference"
##  Oper="AlgebraicExtension" />
##  <P/> 
##  The  underlying  tower of  finite  fields  can  be constructed  with  <Ref
##  Func="StandardFiniteFieldTower" /> or it can  be accessed from <C>F</C> by
##  <Ref Attr="Tower" />.
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
##  REMARK:  Note that  in this  version of  the <Package>StandardFF</Package>
##  package  the  fields returned  by  <Ref  Func="FF"/>  are not  cached  (in
##  contrast  to the  &GAP; function  <Ref BookName="Reference"  Func="GF"/>).
##  Furthermore, currently the package  only supports arithmetic operations of
##  elements in the same field (and  the prime field, for convenience). To add
##  or  multiply elements  in  different  extensions of  the  prime field  use
##  explicit embeddings (see <Ref  Meth="Embedding" Label="for standard finite
##  fields"/>).
##  <Example>gap> Fp := FF(2, 1);
##  GF(2)
##  gap> F := FF(2, 100);
##  FF(2, 100)
##  gap> StandardFiniteFieldTower(2, 100);
##  FFTower(2;2,2,5,5)
##  gap> T := Tower(F);
##  FFTower(2;2,2,5,5)
##  gap> Size(F);
##  1267650600228229401496703205376
##  gap> Size(F) = Size(T);
##  true
##  gap> p := NextPrimeInt(10^50);
##  100000000000000000000000000000000000000000000000151
##  gap> K := FF(p, 60);
##  FF(100000000000000000000000000000000000000000000000151, 60)
##  gap> LogInt(Size(K), 10);
##  3000
##  gap> F := FF(13, 9*5);
##  FF(13, 45)
##  gap> StandardPrimeDegreePolynomial(13, 3, 1);
##  x3_1^3+Z(13)^10
##  gap> StandardPrimeDegreePolynomial(13, 3, 2);
##  x3_2^3-x3_1
##  gap> StandardPrimeDegreePolynomial(13, 5, 1);
##  x5_1^5+Z(13)^3*x5_1-Z(13)^0
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
##  
InstallGlobalFunction(StandardPrimeDegreePolynomial, 
function(p, r, k, variant...)
  if not IsPrimeInt(p) then
    Error("StandardPrimeDegreePolynomial: characteristic p must be a prime.");
  elif not IsPrimeInt(r) then
    Error("StandardPrimeDegreePolynomial: r must be a prime.");
  fi;
  if p = r then
    return SFFHelper.ArtinSchreierStandardFF(p, k)[1];
  elif 2 = r then
    return SFFHelper.2PowerStandardFF(p, k)[1];
  else
    if "Trace" in variant then
      return SFFHelper.rPowerStandardFFTrace(p, r, k)[1];
    else
      return SFFHelper.rPowerStandardFFRandom(p, r, k)[1];
    fi;
  fi;
end);

InstallGlobalFunction(StandardFiniteFieldTower, function(p, n, variant...)
  local pl, f, pols, res, primposs, a, i;
  if n=1 then 
    res := GF(p);
    SetIsStandardPrimeField(res, true);
    return res; 
  fi;
  if "Trace" in variant then
    pl := true;
  else
    pl := false;
  fi;
  f := Collected(FactorsInt(n));
  pols := [];
  for a in f do
    for i in [1..a[2]] do
      if pl then
        Add(pols, StandardPrimeDegreePolynomial(p, a[1], i, "Trace"));
      else
        Add(pols, StandardPrimeDegreePolynomial(p, a[1], i));
      fi;
    od;
  od;
  res := FiniteFieldTower(p, pols);
  primposs := [f[1][2]];
  for i in [2..Length(f)] do
    Add(primposs, primposs[i-1]+f[i][2]);
  od;
  res!.primposs := primposs;
  if not pl then
    SetFilterObj(res, IsStandardFiniteFieldTower);
  fi;
  SetLeftActingDomain(res, FF(p,1));
  return res;
end);

# standardized primitive element is sum of last prime power generators
InstallMethod(PrimitiveElement, ["IsStandardFiniteFieldTower"], function(L)
  local gen;
  gen := GeneratorsOfField(L);
## Would also work:   return Sum(gen{L!.primposs});
  return Product(gen{L!.primposs});
end);

# caching minimal polynomials of primitive elements
# (as Steinitz number of poly minus leading term)
SFFHelper.MinimalPolynomialCACHE := function(Fp, xmat, iv)
  local p, n, res, c;
  if UserPreference("StandardFFUseCache") <> true then
    return MinimalPolynomial(Fp, xmat, iv);
  fi;
  if not IsBound(SFFHelper.MIPOPRIMCACHE) then
    ReadPackage("StandardFF", "data/MIPOPRIMCACHE");
  fi;
  p := Characteristic(Fp);
  n := Length(xmat);
  if p < 10000 then
    if IsBound(SFFHelper.MIPOPRIMCACHE[p]) and 
                   IsBound(SFFHelper.MIPOPRIMCACHE[p][n]) then 
      if IsInt(iv) then
        iv := Indeterminate(Fp, iv);
      fi;
      res := CoefficientsQadic(SFFHelper.MIPOPRIMCACHE[p][n], p);
      while Length(res) < n do
        Add(res, 0);
      od;
      Add(res, 1);
      return ValuePol(One(Fp)*res, iv);
    fi;
  fi;
  res := MinimalPolynomial(Fp, xmat, iv);
  if p < 10000 then
    c := List(CoefficientsOfUnivariatePolynomial(res), IntFFE);
    Remove(c);
    while c[Length(c)] = 0 do
      Remove(c);
    od;
    c := ValuePol(c, p);
    if not IsBound(SFFHelper.MIPOPRIMCACHE[p]) then
      SFFHelper.MIPOPRIMCACHE[p] := [];
    fi;
    SFFHelper.MIPOPRIMCACHE[p][n] := c;
  fi;
  return res;
end;


# standard finite field as simple extension over GF(p)
InstallGlobalFunction(StandardFiniteField, function(p, n)
  local Fp, K, x, vnam, v, iv, xmat, mp, xp, Kp, i, fam;
  if n=1 then 
    Kp := GF(p);
    SetIsStandardPrimeField(Kp, true);
    return Kp; 
  fi;
  Fp := FF(p,1);
  K := StandardFiniteFieldTower(p, n);
  x := PrimitiveElement(K);
  vnam := Concatenation("x", String(n));
  v := Indeterminate(Fp, vnam);
  iv := IndeterminateNumberOfUnivariateLaurentPolynomial(v);
  # x as very sparse matrix acting on tower basis
  xmat := List(TowerBasis(K), b-> AsVector(b*x));
  # minimal polynomial of x is minimal polynomial of xmat
  #mp := MinimalPolynomial(Fp, xmat, iv);
  mp := SFFHelper.MinimalPolynomialCACHE(Fp, xmat, iv);
  # we store the powers of x in tower basis
  xp := [0*xmat[1]];
  xp[1][1] := One(Fp);
  for i in [1..n-1] do
    Add(xp, xp[i] * xmat);
  od;
  Kp := AlgebraicExtensionNC(GF(p), mp, vnam);
  Setter(Tower)(Kp, K);
  Setter(IsStandardFiniteField)(Kp, true);
  Setter(PrimitivePowersInTowerBasis)(Kp, xp);
  # let elements know to be in standard field
  fam := FamilyObj(RootOfDefiningPolynomial(Kp));
  fam!.extType := Subtype(fam!.extType, IsStandardFiniteFieldElement);
  fam!.baseType := Subtype(fam!.baseType, IsStandardFiniteFieldElement);
  SetFilterObj(OneImmutable(Kp), IsStandardFiniteFieldElement);
  SetFilterObj(ZeroImmutable(Kp), IsStandardFiniteFieldElement);
  SetFilterObj(RootOfDefiningPolynomial(Kp), IsStandardFiniteFieldElement);

  return Kp;
end);

# default for other fields
InstallOtherMethod(IsStandardFiniteField, ["IsField"], ReturnFalse);

# maps between elements of standard finite field and its
# standard finite field tower
InstallMethod(TowerBasis, ["IsStandardFiniteField"], function(F)
  return PrimitivePowersInTowerBasis(F)^-1;
end);
# list of degrees of the monomials in tower basis
InstallMethod(TowerBasisMap, ["IsStandardFiniteField"], function(F)
  return StdMon(DegreeOverPrimeField(F))[2];
end);
InstallOtherMethod(TowerBasisMap, ["IsStandardFiniteFieldTower"], function(F)
  return StdMon(DegreeOverPrimeField(F))[2];
end);

InstallMethod(AsVector, ["IsStandardFiniteFieldElement"], function(x)
  local F;
  F := FamilyObj(x)!.wholeField;
  return ExtRepOfObj(x) * PrimitivePowersInTowerBasis(F);
end);

InstallMethod(AsPolynomial, ["IsStandardFiniteFieldElement"], function(x)
  return PolynomialFiniteFieldTowerElement(ToTowerElement(DefaultField(x), x));
end);
InstallOtherMethod(AsPolynomial, ["IsFFE"], IdFunc);

InstallOtherMethod(AsVector, 
    ["IsStandardFiniteField", "IsStandardFiniteFieldElement"], 
function(F, x)
  return ExtRepOfObj(x) * PrimitivePowersInTowerBasis(F);
end);

InstallMethod(ToTowerElement, 
     ["IsStandardFiniteField", "IsStandardFiniteFieldElement"],
function(F, c)
  return FiniteFieldTowerElementVector(Tower(F), AsVector(F, c));
end);

InstallOtherMethod(ToTowerElement, ["IsStandardPrimeField", "IsFFE"],
function(F, c)
  return c;
end);


##  <#GAPDoc Label="IsFF">
##  <ManSection>
##  <Heading>Filters for standard fields</Heading>
##  <Prop Name="IsStandardPrimeField" Arg="F" />
##  <Prop Name="IsStandardFiniteField" Arg="F" />
##  <Filt Name="IsStandardFiniteFieldElement" Arg="x" Type="Category"/>
##  <Filt Name="IsStandardFiniteFieldTower" Arg="T" Type="Category"/>
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
##  And  fields   constructed  by   <C>StandardFiniteFieldTower(p,k)</C>  with
##  <C>k</C> <M>&gt; 1</M> (that is proper extensions of the prime field) have
##  the property <Ref Prop="IsStandardFiniteFieldTower"/>.
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
##  gap> IsStandardFiniteField(Tower(F));
##  false
##  gap> IsStandardFiniteFieldTower(Tower(F));
##  true
##  gap> IsStandardFiniteFieldTower(StandardFiniteFieldTower(11,2));
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
InstallOtherMethod(FromTowerElement, ["IsStandardPrimeField", "IsInt"],
function(F, nr)
  return nr*One(F);
end);

##  <#GAPDoc Label="FFElmConv">
##  <ManSection>
##  <Heading>Maps for elements of standard finite fields</Heading>
##  <Oper Name="ToTowerElement" Arg="F, a"/>
##  <Returns>an element in <C>Tower(F)</C> </Returns>
##  <Oper Name="FromTowerElement" Arg="F, t"/>
##  <Returns>an element in <C>F</C> </Returns>
##  <Meth Name="AsVector" Label="for elements in standard finite fields"
##  Arg="a"/>
##  <Returns>a vector over prime field of <C>F</C> </Returns>
##  <Meth Name="ElementVector" Arg="F, v"/>
##  <Returns>an element in <C>F</C> </Returns>
##  <Meth Name="AsPolynomial" Label="for elements in standard finite fields"
##  Arg="a"/>
##  <Returns>a polynomial in variables of <C>Tower(F)</C> </Returns>
##  <Meth Name="ElementPolynomial" Arg="F, pol"/>
##  <Returns>an element in <C>F</C> </Returns>
##  <Meth Name="SteinitzNumber" Label="for standard finite field elements"
##   Arg="a"/>
##  <Returns>an integer</Returns>
##  <Meth Name="ElementSteinitzNumber" Label="for standard finite fields" 
##  Arg="F, nr"/>
##  <Returns>an element in <C>F</C></Returns>
##  <Description>
##  Here,    <A>F</A>   is    always   a    standard   finite    field   (<Ref
##  Filt="IsStandardFiniteField"/>), <A>a</A>  is an  element of  <A>F</A> and
##  <A>t</A> is an element of <C>Tower(<A>F</A>)</C>.
##  <P/> 
##  Then    <Ref   Oper="ToTowerElement"    />   returns   the   element    in
##  <C>Tower(<A>F</A>)</C>    corresponding    to     <A>a</A>,    and    <Ref
##  Oper="FromTowerElement" /> does the inverse.
##  <P/>
##  <Ref  Meth="AsVector" Label="for  elements  in  standard finite  fields"/>
##  returns the coefficient vector of <A>a</A> with respect to the tower basis
##  of  <C>Tower(<A>F</A>)</C>. And  <Ref  Meth="ElementVector"/> returns  the
##  element of <A>F</A> with the given coefficient vector.
##  <P/>
##  Similarly, <Ref Meth="AsPolynomial" Label="for elements in standard finite
##  fields"/> returns the (reduced)  polynomial in the indeterminates defining
##  <C>Tower(<A>F</A>)</C>  representing  <C>ToTower(<A>F</A>,  <A>a</A>)</C>.
##  And  <Ref  Meth="ElementPolynomial"/>  returns  the  element  of  <A>F</A>
##  represented by the given polynomial (which does not need to be reduced).
##  <P/>
##  Finally,    <Ref    Meth="SteinitzNumber"/>     returns    the    Steinitz
##  number     of     <C>ToTower(<A>F</A>,     <A>a</A>)</C>.     And     <Ref
##  Meth="ElementSteinitzNumber"/>  returns the  element  with given  Steinitz
##  number.
##  <Example><![CDATA[gap> F := FF(17, 12);
##  FF(17, 12)
##  gap> T := Tower(F);
##  FFTower(17;2,2,3)
##  gap> a := PrimitiveElement(F);; a := a^11-3*a^5+a;
##  x12^11+Z(17)^9*x12^5+x12
##  gap> t := ToTowerElement(F, a);
##  <FFTower(17;2,2,3):[ 0, 0, 10, 16, 0, 0, 3, 15, 0, 0, 16, 9 ]>
##  gap> a^12345 = FromTowerElement(F, t^12345);
##  true
##  gap> v := AsVector(a);
##  < mutable compressed vector length 12 over GF(17) >
##  gap> a = ElementVector(F, v);
##  true
##  gap> ExtRepOfObj(a) = v * TowerBasis(F);
##  true
##  gap> pol := AsPolynomial(a);;
##  gap> ElementPolynomial(F, pol^10) = a^10;
##  true
##  gap> nr := SteinitzNumber(a);
##  340709196750181
##  gap> a = ElementSteinitzNumber(F, nr);
##  true
##  gap> rgens := GeneratorsOfField(T); # generators of prime degree steps
##  [ <FFTower(17;2,2,3):[ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]>, 
##    <FFTower(17;2,2,3):[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]>, 
##    <FFTower(17;2,2,3):[ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ]> ]
##  gap> FromTowerElement(F, rgens[2]*rgens[3]) = PrimitiveElement(F);
##  true
##  gap> ## primitive element of FF(17, 6)
##  gap> y := FromTowerElement(F, rgens[1] * rgens[3]);;
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
InstallMethod(FromTowerElement, ["IsStandardFiniteField", "IsRingElement"],
function(F, tc)
  return ElementVector(F, AsVector(tc));
end);
InstallMethod(ElementPolynomial, ["IsStandardFiniteField", "IsPolynomial"],
function(F, pol)
  return FromTowerElement(F, FiniteFieldTowerElement(Tower(F), pol));
end);

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
InstallMethod(ViewString, ["IsStandardFiniteFieldTower"], function(F)
  local res;
  if not IsFiniteFieldTower(F!.prev) then
    return Concatenation("FFTower(", String(Size(F!.prev)),
                                ";",String(F!.deg),")");
  fi;
  res := Concatenation("FF", ViewString(F!.prev));
  Remove(res);
  Append(res,",");
  Append(res, String(F!.deg));
  Append(res,")");
  return res;
end);

##  InstallMethod(PrimitivePowersInTowerBasis, ["IsStandardFiniteField"],
##  function(F)
##    local K, p, x, tdeg, z, xpowers, vec, i;
##    K := Tower(F);
##    p := Characteristic(K);
##    x := PrimitiveElement(K);
##    tdeg := DegreeOverPrimeField(K);
##    z := Zero(GF(p));
##    xpowers := [One(K)];
##    for i in [1..tdeg-1] do
##      Add(xpowers, xpowers[i]*x);
##    od;
##    for i in [1..tdeg] do
##      vec := AsVector(xpowers[i]);
##      while Length(vec) < K!.deg do
##        Add(vec, z);
##      od;
##      xpowers[i] := vec;
##    od;
##    ConvertToMatrixRep(xpowers, p);
##    return xpowers;
##  end);

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
smd := function(n)
  local f, a, res, new, i;
  if n = 1 then return [1]; fi;
  f := Collected(Factors(n));
  a := f[Length(f)];
  res := smd(n/a[1]);
  new := List(res, l-> LcmInt(l, a[1]^a[2]));
  for i in [1..a[1]-1] do
    Append(res, new);
  od;
  return res;
end;

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
##  x45
##  gap> x := y^emb;;
##  gap> ((y+One(H))^12345)^emb = (x+One(F))^12345;
##  true
##  gap> PreImageElm(emb, x^5);
##  x45^5
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
  local p, d, a, monK, monL, map, z, famL, famK, fun, pre, emb;
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
##  The    argument     <A>a</A>    must    be    an     element    in    <Ref
##  Filt="IsStandardFiniteFieldElement"/> or  an element in a  standard finite
##  field tower. Then <Ref Oper="SteinitzPair"/> returns a pair <C>[d, nr]</C>
##  where <C>d</C>  is the degree  of <A>a</A>  over the prime  field <C>FF(p,
##  1)</C>  and <C>nr</C>  is the  Steinitz number  of <A>a</A>  considered as
##  element of <C>FF(p, d)</C>.
##  <P/>
##  In the  second variant a  standard finite field  or its tower  <A>K</A> is
##  given and the Steinitz number of an  element in <A>K</A> and the result is
##  the Steinitz pair of the corresponding element.
##  <P/>
##  The inverse  map is  provided by a  method for  <Ref Meth="SteinitzNumber"
##  Label="for  Steinitz pair"/>  which gets  a standard  finite field  or its
##  tower and a Steinitz pair.
##  <Example>gap> F := FF(7, 360);
##  FF(7, 360)
##  gap> T := Tower(F);
##  FFTower(7;2,2,2,3,3,5)
##  gap> rgen := GeneratorsOfField(T);;
##  gap> t := rgen[2] * rgen[4];; # prim. elt of FF(7,12)
##  gap> sp := SteinitzPair(t);
##  [ 12, 117649 ]
##  gap> a := FromTowerElement(F, t);;
##  gap> H := FF(7, 12);
##  FF(7, 12)
##  gap> ElementSteinitzNumber(H, 117649);
##  x12
##  gap> b := ElementSteinitzNumber(H, 117649);
##  x12
##  gap> Value(MinimalPolynomial(FF(7,1), a), b);
##  !0*Z(7)
##  gap> nr := SteinitzNumber(a);
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
  map := TowerBasisMap(K);
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
InstallOtherMethod(SteinitzPair, 
    ["IsStandardFiniteFieldTower", "IsFiniteFieldTowerElement"], 
    GAPInfo.tmpmeth);
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
InstallOtherMethod(SteinitzPair, 
    ["IsStandardFiniteFieldTower", "IsInt"], GAPInfo.tmpmeth);
Unbind(GAPInfo.tmpmeth);
InstallOtherMethod(SteinitzPair, ["IsStandardPrimeField", "IsInt"], 
function(F, nr)
  return [1, nr];
end);

InstallMethod(SteinitzPair, ["IsStandardFiniteFieldElement"], 
function(x)
  return SteinitzPair(FamilyObj(x)!.wholeField, x);
end);
InstallMethod(SteinitzPair, ["IsFiniteFieldTowerElement"], 
function(x)
  if not IsStandardFiniteFieldTower(x!.field) then
    Error("Use SteinitzPair only for standard finite fields");
    return fail;
  fi;
  return SteinitzPair(x!.field, x);
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
  map := TowerBasisMap(K);
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
InstallOtherMethod(SteinitzNumber, ["IsStandardFiniteFieldTower", "IsList"],
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
##  <Example>gap> F := FF(23,18);
##  FF(23, 18)
##  gap> st := SteinitzPairConwayGenerator(F);
##  [ 18, 2553999235694638106966959 ]
##  gap> st9 := SteinitzPairConwayGenerator(FF(23,9));
##  [ 9, 841141205152 ]
##  gap> st6 := SteinitzPairConwayGenerator(FF(23,6));
##  [ 6, 36792796 ]
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
##  [ 2553999235694638106966959 ]
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
##  The argument <A>F</A> must be  a standard finite field, say <C>FF(p,n)</C>
##  such  that &GAP;  can generate  <C>GF(p,n)</C>. This  function returns  an
##  isomorphism of fields from <C>GF(p,n)</C> to <A>F</A>.
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
