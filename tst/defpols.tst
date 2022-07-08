gap> l := [];; c := [];;
gap> for pd in [[2,[1..32]],[7,[1..30]], [131, [2,3,5,9,11,2*5*11]],
>            [617,[2,3,4,5,6,7,8,9,2*7*5]], [36893488147419103363,
>            [1,2,3,4,5,6,2*3*11]]] do
>   for d in pd[2] do
>     F := FF(pd[1], d);
>     pl := DefiningPolynomial(F);
>     pl := List(CoefficientsOfUnivariatePolynomial(pl), IntFFE);
>     Add(l, pl);
>     #if Log2Int(Size(F)) < 3000 then # disabled due to bugs in gap dev
>     if pd[1] <> 36893488147419103363 then
>       pl := MinimalPolynomial(FF(pd[1],1), StandardPrimitiveRoot(F));
>       pl := List(CoefficientsOfUnivariatePolynomial(pl), IntFFE);
>       Add(c, pl);
>     fi;
>   od;
> od;
gap> ReadPackage("StandardFF", "tst/defpols.res");
true
gap> l = lcompare;
true
gap> c = ccompare{[1..Length(c)]};
true
