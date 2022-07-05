gap> START_TEST("Standard Cyc");
gap> F := FF(2,48);
FF(2, 48)
gap> x := StandardPrimitiveRoot(F);;
gap> for r in Factors(Size(F)-1) do 
> Print(r, ":  ", SteinitzPair(StandardCyclicGenerator(F, r)), " "); 
> Print(x^((Size(F)-1)/r) = StandardCyclicGenerator(F, r), "\n");
> od;
3:  [ 2, 3 ] true
3:  [ 2, 3 ] true
5:  [ 4, 13 ] true
7:  [ 3, 2 ] true
13:  [ 12, 2991 ] true
17:  [ 8, 72 ] true
97:  [ 48, 32548224419439 ] true
241:  [ 24, 15239999 ] true
257:  [ 16, 34287 ] true
673:  [ 48, 22400786130061 ] true
gap> st := SteinitzPairConwayGenerator(F);
[ 48, 80561746637265 ]
gap> lg := DLog(x, ElementSteinitzNumber(F, st[2]));
99786161807144
gap> x^lg = ElementSteinitzNumber(F, st[2]);
true
gap> F := FF(5,12);
FF(5, 12)
gap> x := StandardPrimitiveRoot(F);;
gap> for r in Factors(Size(F)-1) do 
> Print(r, ":  ", SteinitzPair(StandardCyclicGenerator(F, r)), " "); 
> Print(x^((Size(F)-1)/r) = StandardCyclicGenerator(F, r), "\n");
> od;
2:  [ 1, 4 ] true
2:  [ 1, 4 ] true
2:  [ 1, 4 ] true
2:  [ 1, 4 ] true
3:  [ 2, 17 ] true
3:  [ 2, 17 ] true
7:  [ 6, 7977 ] true
13:  [ 4, 467 ] true
31:  [ 3, 49 ] true
601:  [ 12, 11502093 ] true
gap> st := SteinitzPairConwayGenerator(F);
[ 12, 112364098 ]
gap> lg := DLog(x, ElementSteinitzNumber(F, st[2]));
202647317
gap> x^lg = ElementSteinitzNumber(F, st[2]);
true
gap> F := FF(271,10);
FF(271, 10)
gap> x := StandardPrimitiveRoot(F);;
gap> for r in Factors(Size(F)-1) do 
> Print(r, ":  ", SteinitzPair(StandardCyclicGenerator(F, r)), " "); 
> Print(x^((Size(F)-1)/r) = StandardCyclicGenerator(F, r), "\n");
> od;
2:  [ 1, 270 ] true
2:  [ 1, 270 ] true
2:  [ 1, 270 ] true
2:  [ 1, 270 ] true
2:  [ 1, 270 ] true
3:  [ 1, 242 ] true
3:  [ 1, 242 ] true
3:  [ 1, 242 ] true
5:  [ 1, 10 ] true
5:  [ 1, 10 ] true
11:  [ 10, 248950766683767957073349 ] true
17:  [ 2, 26625 ] true
31:  [ 10, 1390116175172035936451651 ] true
41:  [ 10, 808456185040088236353230 ] true
61:  [ 10, 350738960767315597820591 ] true
251:  [ 5, 1348573917257 ] true
6301:  [ 10, 506402625548764727439891 ] true
4313591:  [ 5, 981986828984 ] true
gap> st := SteinitzPairConwayGenerator(F);
[ 10, 537050755118426958133182 ]
gap> lg := DLog(x, ElementSteinitzNumber(F, st[2]));
462284005351058468777461
gap> x^lg = ElementSteinitzNumber(F, st[2]);
true
gap> STOP_TEST("Standard Cyc", 0);

