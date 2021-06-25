gap> START_TEST("FF towers");

# from manual
gap> x := Indeterminate(GF(13), "x");;
gap> y := Indeterminate(GF(13), "y");;
gap> pol1 := x^2 - 2; pol2 := y^2 - x;
x^2+Z(13)^7
y^2-x
gap> # printed by size of base field and degrees of polynomials
gap> F1 := FiniteFieldTower(GF(13), pol1);
Tower(13;2)
gap> PreviousFiniteFieldTower(F1);
GF(13)
gap> F := ExtendFiniteFieldTower(F1, pol2);
Tower(13;2,2)
gap> PreviousFiniteFieldTower(F);
Tower(13;2)
gap> # elements are printed by coefficients wrt. tower basis
gap> gen := GeneratorsOfField(F);
[ <Tower(13;2,2):[ 0, 1, 0, 0 ]>, <Tower(13;2,2):[ 0, 0, 1, 0 ]> ]
gap> Sum(gen); Product(gen);
<Tower(13;2,2):[ 0, 1, 1, 0 ]>
<Tower(13;2,2):[ 0, 0, 0, 1 ]>
gap> Random(F)^(Size(F)-1);
<Tower(13;2,2):[ 1, 0, 0, 0 ]>
gap> Zero(F); One(F);
<Tower(13;2,2):[ 0, 0, 0, 0 ]>
<Tower(13;2,2):[ 1, 0, 0, 0 ]>
gap> # with F as above
gap> TowerBasis(F);
[ <Tower(13;2,2):[ 1, 0, 0, 0 ]>, <Tower(13;2,2):[ 0, 1, 0, 0 ]>, 
  <Tower(13;2,2):[ 0, 0, 1, 0 ]>, <Tower(13;2,2):[ 0, 0, 0, 1 ]> ]
gap> # F as above
gap> pol := x^4 - x^2*y^3 - 2;
-x^2*y^3+x^4+Z(13)^7
gap> c := FiniteFieldTowerElement(F, pol);
<Tower(13;2,2):[ 2, 0, 0, 11 ]>
gap> redpol := PolynomialFiniteFieldTowerElement(c);
Z(13)^7*x*y+Z(13)
gap> c = FiniteFieldTowerElement(F, redpol);
true
gap> # F as above
gap> c := Sum(GeneratorsOfField(F))^123456;
<Tower(13;2,2):[ 6, 7, 0, 10 ]>
gap> vec := AsVector(c);
[ Z(13)^5, Z(13)^11, 0*Z(13), Z(13)^10 ]
gap> c = FiniteFieldTowerElementVector(F, vec);
true
gap> # F as above, base field is GF(13)
gap> List(GeneratorsOfField(F), SteinitzNumber);
[ 13, 169 ]
gap> ElementSteinitzNumber(F, Size(F)-1);
<Tower(13;2,2):[ 12, 12, 12, 12 ]>

# more cases
gap> f := FiniteFieldTower(GF(25));
GF(5^2)
gap> xx := Indeterminate(GF(5), "xx");
xx
gap> f := FiniteFieldTower(GF(5), xx-Z(5));
Tower(5;1)
gap> Size(f);
5
gap> f := FiniteFieldTower(GF(5), xx^3-Z(5)); # not irreducible
Tower(5;3)
gap> a:=ElementSteinitzNumber(f,111);
<Tower(5;3):[ 1, 2, 4 ]>
gap> a^(5^3-1); # nonsense
<Tower(5;3):[ 2, 4, 3 ]>
gap> f := FiniteFieldTower(GF(5), xx^4-Z(5));
Tower(5;4)
gap> a:=ElementSteinitzNumber(f,500);
<Tower(5;4):[ 0, 0, 0, 4 ]>
gap> a^(5^4-1);
<Tower(5;4):[ 1, 0, 0, 0 ]>
gap> for p in [2,7,17] do
> f := GF(p);
> xx := X(f,"xx"); yy := X(f, "yy"); zz := X(f, "zz");
> px := ValuePol(ConwayPol(p,2),xx);
> py := ValuePol(ConwayPol(p,11),yy);
> pz := ValuePol(ConwayPol(p,13),zz);
> F := FiniteFieldTower(f, px);
> F := ExtendFiniteFieldTower(F, py, pz);
> ViewObj([f, BaseFiniteFieldTower(F), F, Size(F) mod 107,
>             PreviousFiniteFieldTower(F)]);
> Print("\n");
> od;
[ GF(2), GF(2), Tower(2;2,11,13), 48, Tower(2;2,11) ]
[ GF(7), GF(7), Tower(7;2,11,13), 4, Tower(7;2,11) ]
[ GF(17), GF(17), Tower(17;2,11,13), 69, Tower(17;2,11) ]
gap> STOP_TEST("FF towers", 0);

