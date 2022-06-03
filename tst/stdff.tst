gap> START_TEST("Standard FF");

# from manual
gap> Fp := FF(2, 1);
GF(2)
gap> F := FF(2, 100);
FF(2, 100)
gap> Size(F);
1267650600228229401496703205376
gap> p := NextPrimeInt(10^50);
100000000000000000000000000000000000000000000000151
gap> K := FF(p, 60);
FF(100000000000000000000000000000000000000000000000151, 60)
gap> LogInt(Size(K), 10);
3000
gap> F := FF(13, 9*5);
FF(13, 45)
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

# mixed arithmetic
gap> F := FF(5,30);
FF(5, 30)
gap> x := StandardPrimitiveRoot(F)^4627462;
ZZ(5,30,[1,2,2,2,0,2,2,1,3,1,4,4,3,2,2,2,4,4,2,1,3,4,3,3,2,1,3,2,0,1])
gap> l := [];;
gap> for k in DivisorsInt(30) do Add(l, x^((5^30-1)/(5^k-1))); od;
gap> for a in l do for b in l do
> erg := [a*b, a+b, a-b, -a+b, a/b, b/a]; od;od;
gap> l[3] * l[5];
ZZ(5,30,[1,4,1,0,3,1,3,0,4,0,0,0,1,0,0,0,3,0,4,0,3,0,4,0,0,4,0,0,2,0])
gap> l1 := List(l, MoveToSmallestStandardField);;
gap> for a in l1 do for b in l do
> erg := [a*b, a+b, a-b, -a+b, a/b, b/a]; od;od;
gap> l1[3] * l[2];
ZZ(5,30,[0,4,0,1,1,1,4,0,2,0,1,3,3,0,3,4,0,3,3,3,3,2,3,3,2,4,4,2,1,0])
gap> for a in l1 do for b in l1 do
> erg := [a*b, a+b, a-b, -a+b, a/b, b/a]; od;od;
gap> for a in l1 do Print(a * l1[4], "\n"); od;
ZZ(5,5,[1,3,1,1,4])
ZZ(5,10,[1,2,4,0,2,3,1,2,1,4])
ZZ(5,15,[4,4,2,4,3,4,2,3,0,0,4,0,1,0,2])
ZZ(5,5,[2,1,3,2,1])
ZZ(5,30,[3,0,2,1,4,2,0,1,3,4,4,0,4,0,4,3,0,1,4,2,1,3,1,3,1,0,2,0,0,2])
ZZ(5,10,[1,0,4,0,0,0,2,3,1,1])
ZZ(5,15,[0,4,0,3,4,1,3,1,0,2,4,1,0,3,3])
ZZ(5,30,[0,3,3,0,2,2,2,3,3,3,4,3,2,2,1,3,2,1,4,2,1,0,1,0,0,2,1,1,0,2])
gap> Z(5)^2 + l[3];
ZZ(5,30,[4,0,0,0,4,0,1,0,3,0,4,0,2,0,2,0,0,0,2,0,2,0,2,0,3,0,1,0,4,0])
gap> Z(5)^2 * l[3];
ZZ(5,30,[0,0,0,0,1,0,4,0,2,0,1,0,3,0,3,0,0,0,3,0,3,0,3,0,2,0,4,0,1,0])
gap> l[6] * Z(5)^2;
ZZ(5,30,[2,2,2,3,2,1,1,3,1,1,4,3,3,4,0,4,3,1,0,0,1,0,3,3,4,3,2,4,1,4])
gap> l[6] + Z(5)^2;
ZZ(5,30,[2,3,3,2,3,4,4,2,4,4,1,2,2,1,0,1,2,4,0,0,4,0,2,2,1,2,3,1,4,1])
gap> STOP_TEST("Standard FF", 0);

