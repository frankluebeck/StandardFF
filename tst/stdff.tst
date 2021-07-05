gap> START_TEST("Standard FF");

# from manual
gap> Fp := FF(2, 1);
GF(2)
gap> F := FF(2, 100);
FF(2, 100)
gap> StandardFiniteFieldTower(2, 100);
FFTower(2;2,2,5,5)
gap> T := Tower(F);
FFTower(2;2,2,5,5)
gap> Size(F);
1267650600228229401496703205376
gap> Size(F) = Size(T);
true
gap> p := NextPrimeInt(10^50);
100000000000000000000000000000000000000000000000151
gap> K := FF(p, 60);
FF(100000000000000000000000000000000000000000000000151, 60)
gap> LogInt(Size(K), 10);
3000
gap> F := FF(13, 9*5);
FF(13, 45)
gap> StandardPrimeDegreePolynomial(13, 3, 1);
x3_1^3+Z(13)^10
gap> StandardPrimeDegreePolynomial(13, 3, 2);
x3_2^3-x3_1
gap> StandardPrimeDegreePolynomial(13, 5, 1);
x5_1^5+Z(13)^3*x5_1-Z(13)^0
gap> F := FF(19,1);
GF(19)
gap> IsStandardFiniteField(F);
false
gap> IsStandardPrimeField(F);
true
gap> F := FF(23,48);
FF(23, 48)
gap> IsStandardFiniteField(F);
true
gap> IsStandardFiniteFieldElement(Random(F));
true
gap> IsStandardFiniteField(Tower(F));
false
gap> IsStandardFiniteFieldTower(Tower(F));
true
gap> IsStandardFiniteFieldTower(StandardFiniteFieldTower(11,2));
true
gap> F := FF(17, 12);
FF(17, 12)
gap> T := Tower(F);
FFTower(17;2,2,3)
gap> a := PrimitiveElement(F);; a := a^11-3*a^5+a;
x12^11+Z(17)^9*x12^5+x12
gap> t := ToTowerElement(F, a);
<FFTower(17;2,2,3):[ 0, 0, 10, 16, 0, 0, 3, 15, 0, 0, 16, 9 ]>
gap> a^12345 = FromTowerElement(F, t^12345);
true
gap> v := AsVector(a);
< mutable compressed vector length 12 over GF(17) >
gap> a = ElementVector(F, v);
true
gap> ExtRepOfObj(a) = v * TowerBasis(F);
true
gap> pol := AsPolynomial(a);;
gap> ElementPolynomial(F, pol^10) = a^10;
true
gap> nr := SteinitzNumber(a);
340709196750181
gap> a = ElementSteinitzNumber(F, nr);
true
gap> rgens := GeneratorsOfField(T); # generators of prime degree steps
[ <FFTower(17;2,2,3):[ 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]>, 
  <FFTower(17;2,2,3):[ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0 ]>, 
  <FFTower(17;2,2,3):[ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 ]> ]
gap> FromTowerElement(F, rgens[2]*rgens[3]) = PrimitiveElement(F);
true
gap> ## primitive element of FF(17, 6)
gap> y := FromTowerElement(F, rgens[1] * rgens[3]);;
gap> x6 := Indeterminate(FF(17,1), "x6");;
gap> MinimalPolynomial(FF(17,1), y, x6) = DefiningPolynomial(FF(17,6));
true
gap> F := FF(7, 360);
FF(7, 360)
gap> T := Tower(F);
FFTower(7;2,2,2,3,3,5)
gap> rgen := GeneratorsOfField(T);;
gap> t := rgen[2] * rgen[4];; # prim. elt of FF(7,12)
gap> SteinitzPair(t);
[ 12, 117649 ]
gap> a := FromTowerElement(F, t);;
gap> H := FF(7, 12);
FF(7, 12)
gap> ElementSteinitzNumber(H, 117649);
x12
gap> b := ElementSteinitzNumber(H, 117649);
x12
gap> Value(MinimalPolynomial(FF(7,1), a), b);
!0*Z(7)
gap> F := FF(7, 360);
FF(7, 360)
gap> H := FF(7, 45);
FF(7, 45)
gap> emb := Embedding(H, F);
MappingByFunction( FF(7, 45), FF(7, 360), function( x ) ... end )
gap> y := PrimitiveElement(H);
x45
gap> x := y^emb;;
gap> ((y+One(H))^12345)^emb = (x+One(F))^12345;
true
gap> PreImageElm(emb, x^5);
x45^5
gap> PreImageElm(emb, PrimitiveElement(F));
fail
gap> SteinitzNumber(y);
13841287201
gap> SteinitzNumber(x) mod 10^50;
72890819326613454654477690085519113574118965817601
gap> SteinitzPair(x);
[ 45, 13841287201 ]

# more examples 
gap> IsStandardPrimeField(GF(47));
true
gap> for p in [2,5,19,1009] do Print([FF(p,1),FF(p,12),FF(p,19)],"\n"); od;
[ GF(2), FF(2, 12), FF(2, 19) ]
[ GF(5), FF(5, 12), FF(5, 19) ]
[ GF(19), FF(19, 12), FF(19, 19) ]
[ GF(1009), FF(1009, 12), FF(1009, 19) ]
gap> F := FF(47, 38);
FF(47, 38)
gap> Order(ElementSteinitzNumber(F, 100));
2208
gap> Size(F);
3465122572092046296464724059395298458186535561070143621795409889
gap> F := FF(2,18);; H := FF(2,6);; emb:=Embedding(H,F);;
gap> ForAll(H, x-> PreImageElm(emb, x^emb) =x);
true
gap> F := FF(3,15);; H := FF(3,5);; emb:=Embedding(H,F);;
gap> ForAll(H, x-> PreImageElm(emb, x^emb) =x);
true
gap> F := FF(17,12);; H := FF(17,1);; emb:=Embedding(H,F);;
gap> ForAll(H, x-> PreImageElm(emb, x^emb) =x);
true
gap> for F in [FF(107,1), FF(7,6), FF(2,13)] do Print(ForAll(F,
> x-> x = ElementSteinitzNumber(F,SteinitzNumber(x))),"\n"); od;
true
true
true
gap> for F in [FF(7,6), FF(2,13)] do T := Tower(F);Print(ForAll(T,
> x-> x = ElementSteinitzNumber(T,SteinitzNumber(x))),"\n"); od;
true
true
gap> STOP_TEST("Standard FF", 0);

