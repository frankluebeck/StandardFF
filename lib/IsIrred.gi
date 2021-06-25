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

# for general polynomials, test for divisors of all lower degrees 
InstallGlobalFunction(IsIrreducibleCoeffList, function(cs, q)
  local d, v, z, o, vq, mat, vv, i;
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
  # precompute in rows of mat: 1^q, X^q, (X^2)^q, ... (x^(d-1))^q mod cs
  mat := [0*cs{[1..d]}];
  mat[1][1] := v[2];
  Add(mat, vq);
  for i in [2..d-1] do
    vv := ProductCoeffs(mat[i], vq);
    ReduceCoeffs(vv, cs);
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
  # now we get further q-powers by matrix multiplication
  # in step i (<=d/2) we check for irreducible factors of degree i
  for i in [2..QuoInt(d,2)] do
    vq := vq*mat;
    while Length(vq) > 0 and IsZero(vq[Length(vq)]) do
      Unbind(vq[Length(vq)]);
    od;
    if Length(GcdCoeffs(vq-v, cs)) > 1 then
      Info(InfoStandardFF, 2, i,"\c");
      return false;
    fi;
  od;
  return true;
end);

# About 1/d of all monic polynomials of degree d with constant term a are
# irreducible.
# We try to find sparse polynomials. After each d tries we allow 
# additional non-zero coefficients.
InstallGlobalFunction(StandardIrreducibleCoeffList, function(K, d, a)
  local felm, l, q, count, dd, seed, rm, k, inc;
  # l is coefficient list, monic and constant coeff a
  l := NullMat(1, d+1, K)[1];
  l[d+1] := One(K);
  l[1] := a;
  l[2] := One(K);
  q := Size(K);
  count := 0;
  # inc is the expected number on non-zero coeffs
  inc := 1;
  while q^inc < 2*d do
    inc := inc+1;
  od;
  dd := 1;
  seed := 1;
  while not IsIrreducibleCoeffList(l, q) do
    # after every d attempts allow inc more non-zero coeffs
    ####if count mod (2*d) = 0 and dd < d then
    if count mod (d) = 0 and dd < d then
      dd := dd + inc;
      if dd > d then
        dd := d;
      fi;
      Info(InfoStandardFF, 2, "(dd=",dd,")\c");
    fi;
    for k in [2..dd] do
      rm := SimpleRandomRange(Size(K), seed);
      seed := rm[2];
      l[k] := ElementSteinitzNumber(K, rm[1]);
    od;
    Info(InfoStandardFF, 2, "*\c");
    count := count+1;
  od;
  return l;
end);

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

StandardIrreducibleCoeffListNTL := function(K, d, a)
  local felm, l, q, count, dd, seed, rm, k, inc;
  # l is coefficient list, monic and constant coeff a
  l := NullMat(1, d+1, K)[1];
  l[d+1] := One(K);
  l[1] := a;
  l[2] := One(K);
  q := Size(K);
##    if not IsPrimeInt(q) then
##      # also for small d, say d < 500?
##      return StandardIrreducibleCoeffList(K, d, a);
##    fi;
  count := 0;
  # inc is the expected number on non-zero coeffs
  inc := 1;
  while q^inc < 2*d do
    inc := inc+1;
  od;
  dd := 1;
  seed := 1;
  while not isirrNTL(l, K) do
    # after every d attempts allow inc more non-zero coeffs
    ####if count mod (2*d) = 0 and dd < d then
    if count mod (d) = 0 and dd < d then
      dd := dd + inc;
      if dd > d then
        dd := d;
      fi;
      Info(InfoStandardFF, 2, "(dd=",dd,")\c");
    fi;
    for k in [2..dd] do
      rm := SimpleRandomRange(Size(K), seed);
      seed := rm[2];
      l[k] := ElementSteinitzNumber(K, rm[1]);
    od;
    Info(InfoStandardFF, 2, "*\c");
    count := count+1;
  od;
  return l;
end;
