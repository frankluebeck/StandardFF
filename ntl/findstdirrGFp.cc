//   g++ -O3 -o findstdirrGFp findstdirrGFp.cc -lntl
#include <NTL/ZZX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/GF2XFactoring.h>

using namespace std;
using namespace NTL;

//  standard affine shift
ZZ sasm(ZZ q)
{
  ZZ m;
  for (m = (4*q)/5; GCD(m, q) > 1; m--);
  return m;
}
ZZ sasa(ZZ q)
{
  return (2*q)/3;
}
ZZ sas(ZZ q, long i, ZZ m, ZZ a)
{
    ZZ res;
    mul(res, m, i);
    add(res, res, a);
    rem(res, res, q);
    return res;
}

int main()
{
   ZZ st, qq, m, a, s;
   long p, inc, d, i, c;
   cin >> p;
   zz_p::init(p);
   long r, count;
   cin >> r;

   for (inc = 1, qq = p; qq < 2*r; qq *= p, inc++);
   d = 0;
   count = 0;
   zz_pX f;
   // first polynomial to test
   SetCoeff(f, 0, -1);
   SetCoeff(f, 1, 1);
   SetCoeff(f, r, 1);
   // the test loop
   for (; ! IterIrredTest(f); count++) {
      //cout << "red: " << count << ":" << f << "\n";
      if ((count % r) == 0) {
        d = d + inc;
        if (d > r-1) d = r-1;
        power(qq, p, d-1);
        m = sasm(qq);
        a = sasa(qq);
        //cout << "d:" << d << "\n";
      }
      s = sas(qq, count, m, a);
      //cout << "st:" << s << "\n";
      for (i=1; i<=d; i++) {
        c = rem(s, p);
        s = (s-c)/p;
        SetCoeff(f, i, c);
      }
   }
   //cout << "# tried " << count << " polynomials.\n";
   //cout << f << "\n";
   SetCoeff(f, r, 0);
   st = 0;
   for (st = 0, i = deg(f); i >= 0; i--) {
     st *= p;
     st += rep(coeff(f, i));
   }
   cout << st << "\n";
}

