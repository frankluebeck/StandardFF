#############################################################################
##
#A  IsIrred.gi                                                   Frank Lübeck
##  
##  The files IsIrred.g{i,d}  contain code to find irreducible
##  polynomials over a finite field, polynomials handled as coefficient
##  lists.
##  

# only for prime degree d
# test that poly divides X^(q^d) - X and is prime to X^q - X
InstallGlobalFunction(IsIrreducibleCoeffListPrimeDegree, function(cs, q)
  local d, v, vp;
  d := Length(cs)-1;
  v := cs{[1,2]};
  v[1] := Zero(v[1]);
  v[2] := One(v[2]);
  vp := PowerModCoeffs(v, q, cs);
  if Length(GcdCoeffs(vp-v, cs)) > 1 then
    return false;
  fi;
  vp := PowerModCoeffs(v, q^d, cs);
  if vp <> v then
    return false;
  fi;
  return true;
end);

# utiltiy for X^(p^k) mod cs - compute (X^i)^p mod cs for i < Length(cs)
# and then use matrix multiplication
InstallGlobalFunction(PMCspecial, function(v, q, cs)
  local p, m, a, vv, pp, i;
  p := SmallestRootInt(q);
  if q = p or 2^Length(cs) > q then
    return PowerModCoeffs(v, q, cs);
  fi;
  m := [[One(v[1])], PowerModCoeffs(v, p, cs)];
  for i in [3..Length(cs)-1] do
    a := ProductCoeffs(m[i-1], m[2]);
    ReduceCoeffs(a, cs);
    ShrinkRowVector(a);
    Add(m, a);
  od;
  vv := v;
  pp := 1;
  while pp < q do
    vv := List(vv, c-> c^p)*m;
    ShrinkRowVector(vv);
    pp := pp*p;
  od;
  return vv;
end);

##  <#GAPDoc Label="IsIrreducibleCoeffList">
##  <ManSection>
##  <Func Name="IsIrreducibleCoeffList" Arg="coeffs, q" />
##  <Returns><K>true</K> or <K>false</K></Returns>
##  <Description>
##  The argument  <A>coeffs</A> must be a  list of elements in  a finite field
##  with <A>q</A> elements (or some subfield of it).
##  <P/>
##  The function checks if the univariate polynomial <M>f</M> with coefficient
##  list <A>coeffs</A>  (ending with  the leading coefficient)  is irreducible
##  over the field with <A>q</A> elements.
##  <P/>
##  The  algorithm  computes the  greatest  common  divisor of  <M>f</M>  with
##  <M>X^{{q^i}} - X</M> for <M>i = 1,  2, \ldots</M> up to half of the degree
##  of <M>f</M>.
##  
##  <Example>gap> cs := Z(3)^0 * ConwayPol(3,8);
##  [ Z(3), Z(3), Z(3), 0*Z(3), Z(3)^0, Z(3), 0*Z(3), 0*Z(3), Z(3)^0 ]
##  gap> IsIrreducibleCoeffList(cs, 3);
##  true
##  gap> F := FF(17,4);; x := PrimitiveElement(F);;
##  gap> cs := [x, x+x^0, 0*x, x^0];
##  [ ZZ(17,4,[0,1,0,0]), ZZ(17,4,[1,1,0,0]), ZZ(17,4,[0]), ZZ(17,4,[1]) ]
##  gap> while not IsIrreducibleCoeffList(cs, 17^4) do
##  >    cs[1] := cs[1] + One(F);
##  > od;
##  gap> cs;
##  [ ZZ(17,4,[8,1,0,0]), ZZ(17,4,[1,1,0,0]), ZZ(17,4,[0]), ZZ(17,4,[1]) ]
##  </Example>
##  </Description>
##  </ManSection>
##  <#/GAPDoc>
# for general polynomials, test for divisors of all lower degrees 
InstallGlobalFunction(IsIrreducibleCoeffList, function(cs, q)
  local d, v, z, o, vq, vqq, mat, vv, m, k, i;
  if not IsMutable(cs) then
    cs := ShallowCopy(cs);
  fi;
  d := Length(cs)-1;
  if IsPrimeInt(d) and q < 5 then
    return IsIrreducibleCoeffListPrimeDegree(cs, q);
  fi;
  v := cs{[1,2]};
  z := Zero(v[1]);
  o := One(v[1]);
  v[1] := z;
  v[2] := o;
 
  # remark: a precomputation of v^i mod cs for i = d..2d-1 and variants
  # of PowerModCoeffs and ReduceCoeffs which use it could be useful
  ##    vq := PowerModCoeffs(v, q, cs);
  vq := PMCspecial(v, q, cs);
  # test for linear factors
  if Length(GcdCoeffs(vq-v, cs)) > 1 then
    Info(InfoStandardFF, 2, 1,"\c");
    return false;
  fi;
  
  # for very small q we take a few more q-th powers,
  # many random polynomials have factors of very small degree
  vqq := vq;
  for k in [1..QuoInt(d, 2*Log2Int(q)+2)] do
    vqq := PowerModCoeffs(vqq, q, cs);
    if Length(GcdCoeffs(vqq-v, cs)) > 1 then
      Info(InfoStandardFF, 2, k+1,"'\c");
      return false;
    fi;
  od;

  # precompute in rows of mat: 1^q, X^q, (X^2)^q, ... (x^(d-1))^q mod cs
  mat := [0*cs{[1..d]}];
  mat[1][1] := v[2];
  Add(mat, vq);
  for i in [2..d-1] do
    vv := ProductCoeffs(mat[i], vq);
    m := ReduceCoeffs(vv, cs);
    if Length(vv) > m then
      vv := vv{[1..m]};
    fi;
    Add(mat, vv);
  od;
  for vv in mat do
    while Length(vv) > d do
      Unbind(vv[Length(vv)]);
    od;
    while Length(vv) < d do
      Add(vv, z);
    od;
  od;

  # Now mat is the Berlekamp matrix, describing the Frobenius
  # on the X^i-basis.
  # We could compute the rank of mat - id (first row is zero),
  # cs is irreducible  iff  the rank is d-1.
  # But when is this good? The loop below stops for random polynomials
  # after few vector-matrix multiplications plus gcd.

  # now we get further q-powers by matrix multiplication
  # in step i (<=d/2) we check for irreducible factors of degree i
##    for i in [2..QuoInt(d,2)] do
  for i in [2..Log2Int(d)] do
    # checking for factors of degree <= Log2(d) first
    # seems a good compromise in experiments
    vq := vq*mat;
    while Length(vq) > 0 and IsZero(vq[Length(vq)]) do
      Unbind(vq[Length(vq)]);
    od;
    if Length(GcdCoeffs(vq-v, cs)) > 1 then
      Info(InfoStandardFF, 3, i,"\c");
      return false;
    fi;
  od;
##    return true;
  # if no small degree factor found we now compute rank
  for i in [1..d] do 
     mat[i][i]:=mat[i][i]-One(mat[i][i]); 
  od;
  return RankMat(mat) = d-1;
end);

# About 1/r of all monic polynomials of degree r with constant term a are
# irreducible.
# We try to find sparse polynomials. After each r tries we allow 
# additional non-zero coefficients.
# See Algorithm 5.5 of the article StdFFCyc.
InstallGlobalFunction(StandardIrreducibleCoeffList, function(K, r, a)
  local l, q, count, inc, d, st, k, qq;
  # l is coefficient list, monic and constant coeff a
  # first polynomial to try is X^r + X + a
  l := NullMat(1, r+1, K)[1];
  l[r+1] := One(K);
  l[1] := a;
  l[2] := One(K);
  q := Size(K);
  # inc is the expected number on non-zero coeffs
  inc := 1;
  while q^inc < 2*r do
    inc := inc+1;
  od;
  # we allow non-zero coeffs up to position d
  d := 0;
  count := 0;
  while not IsIrreducibleCoeffList(l, q) do
    # after every r attempts allow inc more non-zero coeffs
    Info(InfoStandardFF, 2, "reducible ",count,":",List(l, IntFFE),"\n");    
    if count mod r = 0 and d < r-1 then
      d := d + inc;
      if d >= r then
        d := r-1;
      fi;
      qq := q^(d-1);
      Info(InfoStandardFF, 2, "(d=",d,")\c");
    fi;
    st := StandardAffineShift(qq, count);
    st := CoefficientsQadic(st, q);
    while Length(st) < d-1 do
      Add(st, 0);
    od;
    for k in [2..d] do
      l[k] := ElementSteinitzNumber(K, st[k-1]);
    od;
    Info(InfoStandardFF, 3, "*\c");
    count := count+1;
  od;
  Info(InfoStandardFF, 2, "found ",count,":",List(l, IntFFE),"\n");    
  return l;
end);

# we leave this undocumented, seems to be interesting only for
# quite large degrees?
isirrNTL := function(cs, K)
  local d, q, prg, inp, res, p, prcoeffs, c, x;
  d:=DirectoriesPackageLibrary("StandardFF","ntl");
  q := Size(K);
  if IsPrimeInt(q) then
    prg := Filename(d, "isirrGFp");
    if prg = fail then
      return fail;
    fi;
    inp := "";
    Append(inp, String(q));
    Append(inp, "\n[\n");
    Append(inp, JoinStringsWithSeparator(List(cs, a-> String(IntFFE(a))), "\n"));
    Append(inp, "\n]\n");
    res := IO_PipeThrough(prg,[],inp);
    if res = "true\n" then 
      return true;
    else
      return false;
    fi;
  else
    prg := Filename(d, "isirrGFq");
    if prg = fail then
      return fail;
    fi;
    p := SmallestRootInt(q);
    if not Size(LeftActingDomain(K)) = p then
      return fail;
    fi;
    inp := "";
    Append(inp, String(p));
    Add(inp, '\n');
    prcoeffs := function(c)
      local a;
      Append(inp, "[\n");
      for a in c do 
        Append(inp, String(IntFFE(a)));
        Add(inp, '\n');
      od;
      Append(inp, "]\n");
    end;

    c := CoefficientsOfUnivariatePolynomial(DefiningPolynomial(K));
    prcoeffs(c);
    Append(inp,"[\n");
    for x in cs do
      c := ShallowCopy(x![1]);
      if not IsList(c) then
        c := [c];
      fi;
      ShrinkRowVector(c);
      prcoeffs(c);
    od;
    Append(inp,"]\n");
    res := IO_PipeThrough(prg,[],inp);
    if res = "true\n" then 
      return true;
    else
      return false;
    fi;
  fi;
end;

