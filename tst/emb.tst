
gap> START_TEST("Steinitz pairs");
gap> F := FF(3,6);
FF(3, 6)
gap> res := true;;
gap> for x in F do 
> nr := SteinitzNumber(x);
> sp := SteinitzPair(x);
> res := res and
> SteinitzNumber(F,sp) = nr and
> SteinitzPair(F,nr) = sp; 
> od;
gap> res;
true
gap> K := FF(2,240);
FF(2, 240)
gap> L := FF(2,10);
FF(2, 10)
gap> emb := Embedding(L, K);;
gap> ForAll(L, x-> SteinitzPair(x) = SteinitzPair(x^emb));
true
gap> ForAll(L, x-> x = PreImageElm(emb, x^emb));
true
gap> PreImageElm(emb, PrimitiveElement(K));
fail

# same but with prime field as subfield
gap> K := FF(17,24);
FF(17, 24)
gap> L := FF(17,1);
GF(17)
gap> emb := Embedding(L, K);;
gap> ForAll(L, x-> SteinitzPair(x) = SteinitzPair(x^emb));
true
gap> ForAll(L, x-> x = PreImageElm(emb, x^emb));
true
gap> PreImageElm(emb, PrimitiveElement(K));
fail
gap> K := FF(499,1);
GF(499)
gap> L := FF(499,1);
GF(499)
gap> emb := Embedding(L, K);;
gap> ForAll(L, x-> SteinitzPair(x) = SteinitzPair(x^emb));
true
gap> ForAll(L, x-> x = PreImageElm(emb, x^emb));
true
gap> PreImageElm(emb, PrimitiveElement(K));
Z(499)
gap> iso := false;;
gap> for pi in [[5,1],[7,2],[13,5]] do
>   F := FF(pi[1], pi[2]);
>   iso := StandardIsomorphismGF(F);
>   K := GF(pi[1]^pi[2]);
>   t := ForAll(List([1..100], i-> Random(K)), x-> x=PreImageElm(iso,x^iso));
>   x := Random(K);
>   Print(pi, t, (x^iso)^4+x^iso = (x^4+x)^iso,"\n");
> od;
[ 5, 1 ]truetrue
[ 7, 2 ]truetrue
[ 13, 5 ]truetrue
gap> ZZ(2,13,Z(2,13)^10);
ZZ(2,13,[0,0,0,0,0,0,0,0,0,0,1,0,0])
gap> ZZ(3,14,Z(3,14));
ZZ(3,14,[2,1,1,0,2,1,0,0,0,1,1,2,0,2])
gap> ZZ(65537,3,Z(65537,3));
ZZ(65537,3,[5776,2318,12323])
gap> STOP_TEST("Steinitz pairs", 0);

