//   g++ -O3 -o findirr findirr.cc -lntl
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/GF2XFactoring.h>

using namespace std;
using namespace NTL;

int main()
{
   ZZ p;
   cin >> p;
   ZZ_p::init(p);
   long n;
   cin >> n;

   ZZ_pX f;
   GF2X g;

   if (p == 2) {
     BuildIrred(g, n);
     cout << g << "\n";
   } else {
     BuildIrred(f, n);
     cout << f << "\n";
   }
}
